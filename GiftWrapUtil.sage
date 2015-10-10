load("HermesNormal.sage")
load("Facet.sage")
load("2dConvexHull.sage")

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
def RemovePts(Pts, PtsToRemove):
	for Pt in PtsToRemove:
		Pts.remove(Pt)
	return Pts

#-------------------------------------------------------------------------------
def CheckNormalFormDim(HNF):
	Dim = 0
	for i in xrange(len(HNF[0])):
		if HNF[Dim][i] != 0:
			Dim += 1
			if Dim == len(HNF):
				break
	return Dim

"""
There are two different ways that we use the hermite normal form. We use it to
get the UCT, which is obtained by stacking points vertically and then calling
echelon_form(). However, we also need to do a more typical HNF in the linear
algebra sense, when we stack the points horizontally. This is the one that is
necessary for getting the normal. For it, include_zero_rows must be false,
because this deals with having an initial set of points that isn't linearly
independent.
"""
#-------------------------------------------------------------------------------
def GetUCT(Pts):
#Stack points vertically
	return (((matrix(matrix(Pts),ZZ).transpose()).echelon_form(include_zero_rows = True, transformation = True)[1])^-1).transpose()

#-------------------------------------------------------------------------------
def GetUCTAndNormal(Vector):
#Stack points vertically
	HNF, UCT = (matrix(matrix(Vector),ZZ).transpose()).echelon_form(include_zero_rows = True, transformation = True)
	#We need to check to make sure we're getting the right
	UCT = ((UCT)^-1)
	if ConvertMatrixToList(UCT.transpose())[0] != Vector:
		print "Unexpected first row after UCT. Should be vector"
		print "Vector = ", Vector
		print "UCT", UCT
		raw_input()
	return UCT.transpose(), ConvertMatrixToList(HNF.transpose())

#-------------------------------------------------------------------------------
def GetHNF(Pts):
	# Instead of thinking about these as points, we should think about these as
	# being vectors.
	if len(Pts) < 2:
		return "Insufficient points for a Hermite Normal Form"

	Vectors = []
	Pt = Pts[0]
	for i in xrange(1,len(Pts)):
		Vectors.append([Pt[j] - Pts[i][j] for j in xrange(len(Pt))])

	HNFMatrix = (matrix(Vectors)).echelon_form(include_zero_rows = False)
	return ConvertMatrixToList(HNFMatrix)

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
def GetHNFFromVectorAndPts(Vector, Pts):
	# Instead of thinking about these as points, we should think about these as
	# being vectors.
	if len(Pts) < 2:
		return "Insufficient points for a Hermite Normal Form"

	Vectors = []
	Pt = Pts[0]
	for i in xrange(1,len(Pts)):
		Vectors.append([Pt[j] - Pts[i][j] for j in xrange(len(Pt))])
	Vectors.append(Vector)
	HNFMatrix = (matrix(Vectors)).echelon_form(include_zero_rows = False)
	return ConvertMatrixToList(HNFMatrix)

#-------------------------------------------------------------------------------
def GetHNFFromVectors(Vectors):
	# Instead of thinking about these as points, we should think about these as
	# being vectors.
	HNFMatrix = (matrix(Vectors)).echelon_form(include_zero_rows = False)
	return ConvertMatrixToList(HNFMatrix)

#-------------------------------------------------------------------------------
def GetNormalFromNonIntVectors(Vectors):
	return list(Matrix(Vectors).echelon_form().transpose().kernel().gens()[0])


#-------------------------------------------------------------------------------
def GetNormalFromHNF(HNF):
	M = matrix(HNF)
	HNF = []
	for Row in M:
		HNF.append(Row)
	Normal = FindNormal(HNF)
	for i in xrange(len(Normal)):
		Normal[i] = Integer(Normal[i])
	return Normal

#-------------------------------------------------------------------------------
def DotProduct(U,V):
	DotProduct = 0
	for i in range(len(U)):
		DotProduct +=  U[i]*V[i]
	return DotProduct

#-------------------------------------------------------------------------------
def TransformPts(Pts, UCT):
	if len(Pts) < 1:
		return []
	TransformMatrix = UCT*(matrix(Pts).transpose())
	return ConvertMatrixToList(TransformMatrix.transpose())

#-------------------------------------------------------------------------------
def TransformPt(Pt, UCT):
	return TransformPts([Pt], UCT)[0]

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
	Vector = []
	for i in xrange(len(Pt1)):
		Vector.append(Pt2[i] - Pt1[i])
	return Vector

#-------------------------------------------------------------------------------
def NormalPointsTowardsPt(Normal, Barycenter, Pt):
	def DistanceBetweenPts(Pt1, Pt2):
		D = 0
		for i in xrange(len(Pt1)):
			D += (Pt1[i] - Pt2[i])^2
		return D
	NegNormal = [-Normal[i] for i in xrange(len(Normal))]
	NormalDistance = DistanceBetweenPts(Barycenter, [Normal[i] + Pt[i] for i in xrange(len(Pt))])
	NegNormalDistance = DistanceBetweenPts(Barycenter, [NegNormal[i] + Pt[i] for i in xrange(len(Pt))])
	if NegNormalDistance > NormalDistance:
		return True
	elif NegNormalDistance < NormalDistance:
		return False
	return

#-------------------------------------------------------------------------------
def FindNewFacetPts(Pts, EdgePts, Normal, KnownFacetPts, NormalThroughFacet):
	MaxAngle = 'Test'
	NewFacetPts = []
	for Pt in Pts:
		if Pt not in KnownFacetPts:
			Vector = MakeVector(EdgePts[0], Pt)
			Denom = DotProduct(Vector, Normal)
			if n(Denom,1000) == 0:
				print "Internal error in FindNewFacetPts, denominator = 0."
				raw_input()
			Num = DotProduct(Vector, NormalThroughFacet)
			Angle = - Num/Denom
			if MaxAngle == 'Test':
				MaxAngle = Angle
				NewFacetPts = [Pt]
				CorrectNum = Num
				CorrectDenom = Denom
			elif Angle < MaxAngle:
				MaxAngle = Angle
				NewFacetPts = [Pt]
				CorrectNum = Num
				CorrectDenom = Denom
			elif Angle == MaxAngle:
				NewFacetPts.append(Pt)
	return NewFacetPts, CorrectNum, CorrectDenom

#-------------------------------------------------------------------------------
def FindNewFacetPtsFromEdge(Pts, EdgePts, Normal, KnownFacetPts):
	# Note that the input normal is an inner normal
	Normal = [Normal[i] for i in xrange(len(Normal))]
	NormalThroughFacet = GetNormalFromHNF(GetHNFFromVectorAndPts(Normal, EdgePts))
	# I want this to be an outer facing normal
	AdditionalPtOnFacet = 'Sentinel'
	for Pt in KnownFacetPts:
		if Pt not in EdgePts:
			AdditionalPtOnFacet = Pt
			break
	# For the sake of consistency, point the normal outwards. Decide later if that's what we want
	if NormalPointsTowardsPt(NormalThroughFacet, AdditionalPtOnFacet, EdgePts[0]):
			NormalThroughFacet = [-NormalThroughFacet[i] for i in xrange(len(NormalThroughFacet))]
	return FindNewFacetPts(Pts, EdgePts, Normal, KnownFacetPts, NormalThroughFacet)

#-------------------------------------------------------------------------------
def CheckAllPtsLieOnOthersideOfFacetHyperplane(Pts, FacetPts, UCT):
	Pts = TransformPts(Pts, UCT)
	for Pt in FacetPts:
		Pts.remove(Pt)
	
	# Start by checking that the UCT sends all of the points to the same side
	VertHyperplaneX = FacetPts[0][0]
	for i in xrange(1,len(FacetPts)):
		if FacetPts[i][0] != VertHyperplaneX:
			print "Internal error, the facet points do not all lie on the same vertical hyperplane after a UCT."
			print FacetPts
			raw_input()

	VertHyperplaneIsGreatest = "Maybe"
	for Pt in Pts:
		if VertHyperplaneIsGreatest == True and Pt[0] > VertHyperplaneX:
			print "Internal error: all points do not lie on the same side of the facet hyperplane"
			print Pt
			print TransformPt(Pt, matrix(UCT^-1, ZZ))
			print VertHyperplaneX
			print Pts
			raw_input()
		elif VertHyperplaneIsGreatest == False and Pt[0] < VertHyperplaneX:
			print "Internal error: all points do not lie on the same side of the facet hyperplane"
			print Pt
			print TransformPt(Pt, matrix(UCT^-1, ZZ))
			print VertHyperplaneX
			print Pts
			raw_input()
		elif VertHyperplaneIsGreatest == "Maybe":
			if Pt[0] < VertHyperplaneX:
				VertHyperplaneIsGreatest = True
			if Pt[0] > VertHyperplaneX:
				VertHyperplaneIsGreatest = False
		elif Pt[0] == VertHyperplaneX:
			print "Internal error: all points do not lie on the same side of the facet hyperplane"
			print Pt
			print TransformPt(Pt, matrix(UCT^-1, ZZ))
			print VertHyperplaneX
			print Pts
			raw_input()
	return

#-------------------------------------------------------------------------------
def WrapMaps(Pts):
	Dim = CheckNormalFormDim(GetHNF(Pts))
	ExistingSPTLPM = "A"
	ExistingLPTSPM = "A"
	Length = len(Pts[0])
	while True:
		Pts, ShortPointToLongPointMap, LongPointToShortPointMap, Dim = MakePointMap(Pts)
		if len(LongPointToShortPointMap.values()[0]) == Length:
			if ExistingSPTLPM == "A":
				ExistingSPTLPM = ShortPointToLongPointMap
				ExistingLPTSPM = LongPointToShortPointMap
			break
		else:
			Length = len(LongPointToShortPointMap.values()[0])
		if ExistingSPTLPM == "A":
			ExistingSPTLPM = ShortPointToLongPointMap
			ExistingLPTSPM = LongPointToShortPointMap
		else:
				ExistingSPTLPM = ComposeMaps(ShortPointToLongPointMap, ExistingSPTLPM)
				ExistingLPTSPM = ComposeMaps(ExistingLPTSPM, LongPointToShortPointMap)
	return Pts, ExistingSPTLPM, ExistingLPTSPM, Dim

#-------------------------------------------------------------------------------
def ComposeMaps(Map1, Map2):
	return {tuple(k): Map2.get(tuple(v)) for k, v in Map1.items()}

#-------------------------------------------------------------------------------
def MakePointMap(Pts):
	HNF = GetHNF(Pts)
	Dim = CheckNormalFormDim(HNF)
	if Dim != len(Pts[0]):
		UCT = GetUCT(GetNormalFromHNF(GetHNF(Pts)))
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
			ShortPointToLongPointMap[tuple(Pt)] = TransformPt(LongPt, ReverseUCT)
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
def PutFacetsInCorrectDimension(Facets, PointMap, Barycenter):
	for Facet in Facets:
		CorrectVertices = []
		for Vertex in Facet.Vertices:
			CorrectVertices.append(PointMap[tuple(Vertex)])
		Facet.Vertices = CorrectVertices
		Normal = GetNormalFromHNF(GetHNF(CorrectVertices))
		if not NormalPointsTowardsPt(Normal, Barycenter, CorrectVertices[0]):
			Normal = [-Normal[i] for i in xrange(len(Normal))]
		Facet.InnerNormal = Normal
	return Facets
	
#-------------------------------------------------------------------------------
def MakeExactUnitVector(Vector):
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

#-------------------------------------------------------------------------------
def CheckFacetAgainstSage(Pts, PtsInFacet):
	Pts.sort()
	PtsInFacet.sort()
	Polytope = Polyhedron(vertices = Pts)
	Facet = Polyhedron(vertices = PtsInFacet)
	FacetList = list(Facet.face_lattice()[-1].ambient_Vrepresentation())
	for i in xrange(len(Polytope.face_lattice())):
		TestFace = list(Polytope.face_lattice()[i].ambient_Vrepresentation())
		TestFace.sort()
		if TestFace == FacetList:
			return
	print "The following points do not lie on a facet: ", PtsInFacet
	raw_input()
	return

#-------------------------------------------------------------------------------
def FindInitialFacet(Pts, Barycenter):
	# Not sure if I ever need the Barycenter. requires some more testing...

	def ScaleNumAndDenom(Num, Denom):
		C = Num^2+Denom^2
		return Num/sqrt(C), Denom/sqrt(C)

	def FindFirstPts(Pts):
		# Just sorting the points
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
	Axes = []
	for i in xrange(2,len(Pts[0])):
		Axes.append([1 if j == i else 0 for j in xrange(len(Pts[0]))])
	while True:
		Dim = CheckNormalFormDim(GetHNF(FirstPts))
		if Dim == len(Pts[0]) - 1:
			break
		ListForNormalThroughFacet = [Normal]
		for i in xrange(1, len(FirstPts)):
			ListForNormalThroughFacet.append(MakeVector(FirstPts[0], FirstPts[i]))
		if len(FirstPts) == 1:
			Dim = 0
		for i in xrange(Dim, len(Axes)):
			ListForNormalThroughFacet.append(Axes[i])
		NormalThroughFacet = MakeExactUnitVector(GetNormalFromNonIntVectors(ListForNormalThroughFacet))
		NewPts, Num, Denom = FindNewFacetPts(Pts, FirstPts, Normal, FirstPts, NormalThroughFacet)
		Num, Denom = ScaleNumAndDenom(Num, Denom)
		FirstPts = NewPts + FirstPts
		PreviousNormal = copy(Normal)
		Normal = [-NormalThroughFacet[i]*Denom + Normal[i]*Num for i in xrange(len(Normal))]
		if not NormalPointsTowardsPt(Normal, Barycenter, FirstPts[0]):
			print "FLIP NORMAL"
			Normal = [-Normal[i] for i in xrange(len(Normal))]
	CheckFacetAgainstSage(Pts, FirstPts)
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
