load("Face.sage")

def PtsAreValid(Pts):
	Dim = len(Pts[0])
	for Pt in Pts:
		if len(Pt) != Dim:
			return False
		for Coord in Pt:
			if type(Coord) is not sage.rings.integer.Integer:
				return False
	return True

#-------------------------------------------------------------------------------
def RemoveDups(Pts):
	UniquePts = []
	for Pt in Pts:
		if Pt not in UniquePts:
			UniquePts.append(Pt)
	return UniquePts

#-------------------------------------------------------------------------------
def CheckNormalFormDim(HNF):
	Dim = 0
	for i in xrange(len(HNF[0])):
		if HNF[Dim][i] != 0:
			Dim += 1
			if Dim == len(HNF):
				break
	return Dim

#-------------------------------------------------------------------------------
def GetHNF(Pts):
	# Instead of thinking about these as points, we should think about these as
	# being vectors.
	if len(Pts) < 2:
		return "Insufficient points for a Hermite Normal Form"

	Vectors = []
	for i in xrange(1,len(Pts)):
		Vectors.append(MakeVector(Pts[0], Pts[i]))

	HNFMatrix = (matrix(Vectors)).echelon_form(include_zero_rows = False)
	return ConvertMatrixToList(HNFMatrix)

#-------------------------------------------------------------------------------
def GetNormalFromHNF(HNF):
	# It feels wrong taking the first generator of the kernel
	return list(Matrix(HNF).transpose().kernel().gens()[0])

#-------------------------------------------------------------------------------
def ConvertMatrixToList(Matrix):
	Pts = []
	for Row in Matrix:
		Pts.append(list(Row))
	for Pt in Pts:
		for i in xrange(len(Pt)):
			Pt[i] = Integer(Pt[i])
	return Pts

#-------------------------------------------------------------------------------
def FindBarycenter(Pts):
	Barycenter = [0 for i in xrange(len(Pts[0]))]
	for i in xrange(len(Pts)):
		for j in xrange(len(Pts[i])):
			Barycenter[j] += Pts[i][j]
	for i in xrange(len(Barycenter)):
		Barycenter[i] = Barycenter[i]/len(Pts)
	return Barycenter

#-------------------------------------------------------------------------------
def MakeVector(Pt1, Pt2):
	# The first point is where the vector originates and the second is where it points
	return [Pt2[i] - Pt1[i] for i in xrange(len(Pt1))]

#-------------------------------------------------------------------------------
def MakeNormalPointInDirectionOfPt(Normal, Barycenter, Pt):
	def DistanceBetweenPts(Pt1, Pt2):
		D = 0
		for i in xrange(len(Pt1)):
			D += (Pt1[i] - Pt2[i])^2
		return D
	NegNormal = [-Normal[i] for i in xrange(len(Normal))]
	NormalDistance = DistanceBetweenPts(Barycenter, [Normal[i] + Pt[i] for i in xrange(len(Pt))])
	NegNormalDistance = DistanceBetweenPts(Barycenter, [NegNormal[i] + Pt[i] for i in xrange(len(Pt))])
	if NegNormalDistance < NormalDistance:
		Normal = [-Normal[i] for i in xrange(len(Normal))]
	return Normal

#-------------------------------------------------------------------------------
def FindNewFacePts(Pts, EdgePts, Normal, KnownFacePts, NormalThroughFace):
	def DotProduct(U,V):
		DotProduct = 0
		for i in range(len(U)):
			DotProduct +=  U[i]*V[i]
		return DotProduct

	MaxAngle = 'Test'
	NewFacePts = []
	for Pt in Pts:
		if Pt not in KnownFacePts:
			Vector = MakeVector(EdgePts[0], Pt)
			Denom = DotProduct(Vector, Normal)
			if n(Denom,1000) == 0:
				print "Internal error in FindNewFacePts, denominator = 0."
				raw_input()
			Num = DotProduct(Vector, NormalThroughFace)
			Angle = - Num/Denom
			if MaxAngle == 'Test':
				MaxAngle = Angle
				NewFacePts = [Pt]
				CorrectNum = Num
				CorrectDenom = Denom
			elif Angle > MaxAngle:
				MaxAngle = Angle
				NewFacePts = [Pt]
				CorrectNum = Num
				CorrectDenom = Denom
			elif Angle == MaxAngle:
				NewFacePts.append(Pt)
	return NewFacePts, CorrectNum, CorrectDenom

#-------------------------------------------------------------------------------
def FindNewFacePtsFromEdge(Pts, EdgePts, Normal, KnownFacePts):
	if len(EdgePts) < 2:
		print "Insufficient points to compute the required Hermite Normal Form"
		raw_input()
	Vectors = [Normal] + [MakeVector(EdgePts[i], EdgePts[0]) for i in xrange(1, len(EdgePts))]
	HNFMatrix = (matrix(Vectors)).echelon_form(include_zero_rows = False)
	NormalThroughFace = GetNormalFromHNF(ConvertMatrixToList(HNFMatrix))
	for Pt in KnownFacePts:
		if Pt not in EdgePts:
			AdditionalPtOnFace = Pt
			break
	# Make the normal an inner normal
	NormalThroughFace = MakeNormalPointInDirectionOfPt(NormalThroughFace, AdditionalPtOnFace, EdgePts[0])
	return FindNewFacePts(Pts, EdgePts, Normal, KnownFacePts, NormalThroughFace)

#-------------------------------------------------------------------------------
def MakeMaps(Pts, ShortPointToLongPointMap = "Start", LongPointToShortPointMap = "Start"):
	def ComposeMaps(Map1, Map2):
		return {tuple(k): Map2.get(tuple(v)) for k, v in Map1.items()}

	while True:
		Length = len(Pts[0])

		Pts, TempShortPointToLongPointMap, TempLongPointToShortPointMap, Dim = MakePointMap(Pts)

		if len(Pts[0]) == Length:
			if ShortPointToLongPointMap == "Start":
				return Pts, TempShortPointToLongPointMap, TempLongPointToShortPointMap, Dim
			else:
				break

		if ShortPointToLongPointMap == "Start":
			ShortPointToLongPointMap = TempShortPointToLongPointMap
			LongPointToShortPointMap = TempLongPointToShortPointMap
		else:
			ShortPointToLongPointMap = ComposeMaps(TempShortPointToLongPointMap, ShortPointToLongPointMap)
			LongPointToShortPointMap = ComposeMaps(LongPointToShortPointMap, TempLongPointToShortPointMap)
			KeysToRemove = []
			for Key in LongPointToShortPointMap:
				if LongPointToShortPointMap[Key] == None:
					KeysToRemove.append(Key)
			for Key in KeysToRemove:
				LongPointToShortPointMap.pop(Key)
	return Pts, ShortPointToLongPointMap, 	LongPointToShortPointMap, Dim

#-------------------------------------------------------------------------------
def MakePointMap(Pts):
	def GetUCT(Pts):
		return (((matrix(matrix(Pts),ZZ).transpose()).echelon_form(include_zero_rows = True, transformation = True)[1])^-1).transpose()

	def TransformPts(Pts, UCT):
		if len(Pts) < 1:
			return []
		TransformMatrix = UCT*(matrix(Pts).transpose())
		return ConvertMatrixToList(TransformMatrix.transpose())

	HNF = GetHNF(Pts)
	Dim = CheckNormalFormDim(HNF)
	if Dim != len(Pts[0]):
		UCT = GetUCT(GetNormalFromHNF(HNF))
		NewPts = TransformPts(Pts, UCT)
		Row = copy(NewPts[0])
		for i in xrange(1, len(NewPts)):
			for j in xrange(len(NewPts[0])):
				if NewPts[i][j] != Row[j]:
					Row[j] = "Diff"
		GoodPts = []
		for i in xrange(len(NewPts)):
			GoodPts.append([])
			for j in xrange(len(NewPts[i])):
				if Row[j] == "Diff":
					GoodPts[i].append(NewPts[i][j])
		ReverseUCT = matrix(UCT^-1, ZZ)
		ShortPointToLongPointMap = {}
		LongPointToShortPointMap = {}
		for Pt in GoodPts:
			LongPt = []
			MyIndex = 0
			for Coord in Row:
				if Coord != "Diff":
					LongPt.append(Coord)
				else:
					LongPt.append(Pt[MyIndex])
					MyIndex += 1
			ShortPointToLongPointMap[tuple(Pt)] = TransformPts(LongPt, ReverseUCT)[0]
			LongPointToShortPointMap[tuple(ShortPointToLongPointMap[tuple(Pt)])] = Pt
	else:
		GoodPts = copy(Pts)
		ShortPointToLongPointMap = {}
		LongPointToShortPointMap = {}
		for Pt in GoodPts:
			ShortPointToLongPointMap[tuple(Pt)] = Pt
			LongPointToShortPointMap[tuple(Pt)] = Pt
	return GoodPts, ShortPointToLongPointMap, LongPointToShortPointMap, Dim

#-------------------------------------------------------------------------------
def FindInitialFacet(Pts):
	def GetUnitNormalThroughFacet(Vectors):
		Vector = list(Matrix(Vectors).echelon_form().transpose().kernel().gens()[0])
		UnitVector = []
		Sum = 0
		for i in Vector:
			Sum += i^2
		ScaleFactor = sqrt(Sum)
		if ScaleFactor == 0 or ScaleFactor == 1:
			return Vector
		elif ScaleFactor == -1:
			return [-Vector[i] for i in xrange(len(Vector))]
		return [Vector[i]/ScaleFactor for i in xrange(len(Vector))]

	def ScaleNumAndDenom(Num, Denom):
		C = Num^2+Denom^2
		return Num/sqrt(C), Denom/sqrt(C)

	def FindFirstPts(Pts):
		FirstPts = []
		FirstPtValue = Pts[0][0]
		for Pt in Pts:
			if Pt[0] == FirstPtValue:
				FirstPts.append(Pt)
			if Pt[0] > FirstPtValue:
				FirstPts = [Pt]
				FirstPtValue = Pt[0]
		return FirstPts

	FirstPts = FindFirstPts(Pts)
	Normal = [-1 if i == 0 else 0 for i in xrange(len(Pts[0]))]
	Axes = [[1 if j == i else 0 for j in xrange(len(Pts[0]))] for i in xrange(2,len(Pts[0]))]
	while True:
		if len(FirstPts) == 1:
			Dim = 0
		else:
			Dim = CheckNormalFormDim(GetHNF(FirstPts))
		if Dim == len(Pts[0]) - 1:
			break
		ListForComputingNormalThroughFacet = [Normal]
		for i in xrange(1, len(FirstPts)):
			ListForComputingNormalThroughFacet.append(MakeVector(FirstPts[0], FirstPts[i]))
		for i in xrange(Dim, len(Axes)):
			ListForComputingNormalThroughFacet.append(Axes[i])
		NormalThroughFacet = GetUnitNormalThroughFacet(ListForComputingNormalThroughFacet)
		NewPts, Num, Denom = FindNewFacePts(Pts, FirstPts, Normal, FirstPts, NormalThroughFacet)
		Num, Denom = ScaleNumAndDenom(Num, Denom)
		FirstPts = NewPts + FirstPts
		PreviousNormal = copy(Normal)
		Normal = [-NormalThroughFacet[i]*Denom + Normal[i]*Num for i in xrange(len(Normal))]
	return FirstPts

#-------------------------------------------------------------------------------
def MakeIndexMaps(Pts):
	PointToIndexMap = {}
	IndexToPointMap = {}
	for i in xrange(len(Pts)):
		Pt = tuple(Pts[i])
		PointToIndexMap[Pt] = i
		IndexToPointMap[i] = Pts[i]
	return PointToIndexMap, IndexToPointMap

#-------------------------------------------------------------------------------
def PrintFaces():
	print "BEGIN"
	for i in xrange(len(Faces)):
		print ""
		print ""
		print ""
		print "----------------------------------------------"
		print "Below are faces of dim", i + 1
		for Face in Faces[i]:
			Face.PrintProps()
	print "END"
	return

#-------------------------------------------------------------------------------
def PrintFaceLens():
	for i in xrange(len(Faces)):
		print "LenFaces[",i,"] = ", len(Faces[i])
	return

#-------------------------------------------------------------------------------
def ConvexHull2d(Pts):
	def FindNextPt(Pts, Pt1):
		def DoTurn(Pt1, Pt2, Pt3):
		 	return cmp((Pt2[0] - Pt1[0])*(Pt3[1] - Pt1[1]) - (Pt3[0] - Pt1[0])*(Pt2[1] - Pt1[1]), 0)
		def Distance(p, q):
			return (q[0] - p[0])^2 + (q[1] - p[1])^2
		Pt2 = Pt1
		for Pt3 in Pts:
			Turn = DoTurn(Pt1, Pt2, Pt3)
			if Turn == -1 or Turn == 0 and Distance(Pt1, Pt3) > Distance(Pt1, Pt2):
				Pt2 = Pt3
		return Pt2
	Hull = [min(Pts)]
	for Pt in Hull:
		Pt2 = FindNextPt(Pts, Pt)
		if Pt2 != Hull[0]:
			Hull.append(Pt2)
	PtsToRemove = []
	for Pt in Pts:
		if Pt not in Hull:
			PtsToRemove.append(Pt)
	return Hull, PtsToRemove

#-------------------------------------------------------------------------------
def CreateCyclicLists(n):
	system = []
	for i in xrange(n-1):
		equation = []
		for j in xrange(n):
			mon = [0 for x in xrange(n)]
			for k in xrange(i+1):
				mon[(j+k)%n]=1
			equation.append(mon)
		system.append(equation)
	mon1 = [0 for x in xrange(n)]
	mon2 = [1 for x in xrange(n)]
	system.append([mon1,mon2])
	return system

#-------------------------------------------------------------------------------
def RemovePts(PtIndicesToRemove, Pts, ShortPointToLongPointMap, LongPointToShortPointMap, IndexToPointMap):
	for PtIndex in PtIndicesToRemove:
		Pt = IndexToPointMap[PtIndex]
		ShortPt = LongPointToShortPointMap.pop(tuple(Pt))
		Pts.remove(ShortPt)
		ShortPointToLongPointMap.pop(tuple(ShortPt))
	return Pts, ShortPointToLongPointMap, LongPointToShortPointMap
