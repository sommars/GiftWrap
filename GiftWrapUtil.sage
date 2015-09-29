load("HermesNormal.sage")
load("2dConvexHull.sage")
load("Facet.sage")

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
	def ConvertMatrixToList(Matrix):
		Pts = []
		for Row in Matrix.transpose():
			Pts.append(list(Row))
		return Pts
#Stack points vertically
	HNF, UCT = (matrix(matrix(Vector),ZZ).transpose()).echelon_form(include_zero_rows = True, transformation = True)
	#We need to check to make sure we're getting the right
	UCT = ((UCT)^-1)
	if ConvertMatrixToList(UCT)[0] != Vector:
		print "Unexpected first row after UCT. Should be vector"
		print "Vector = ", Vector
		print "UCT", UCT
		raw_input()
	return UCT.transpose(), ConvertMatrixToList(HNF)

#-------------------------------------------------------------------------------
def GetHNF(Pts):
	# Instead of thinking about these as points, we should think about these as
	# being vectors.
	def ConvertMatrixToList(Matrix):
		Pts = []
		for Row in Matrix:
			Pts.append(list(Row))
		return Pts
	if len(Pts) < 2:
		return "Insufficient points for a Hermite Normal Form"

	Vectors = []
	Pt = Pts[0]
	for i in xrange(1,len(Pts)):
		Vectors.append([Pt[j] - Pts[i][j] for j in xrange(len(Pt))])

	HNFMatrix = (matrix(Vectors)).echelon_form(include_zero_rows = False)
	return ConvertMatrixToList(HNFMatrix)

#-------------------------------------------------------------------------------
def GetHNFFromVectorAndPts(Vector, Pts):
	# Instead of thinking about these as points, we should think about these as
	# being vectors.
	def ConvertMatrixToList(Matrix):
		Pts = []
		for Row in Matrix:
			Pts.append(list(Row))
		return Pts
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
	def ConvertMatrixToList(Matrix):
		Pts = []
		for Row in Matrix:
			Pts.append(list(Row))
		return Pts
	HNFMatrix = (matrix(Vectors)).echelon_form(include_zero_rows = False)
	return ConvertMatrixToList(HNFMatrix)

#-------------------------------------------------------------------------------
def GetNormalFromNonIntVectors(Vectors):
	LE = len(Matrix(Vectors).echelon_form().transpose().kernel().gens())
	#if LE > 1:
	#	L1 = Matrix(Vectors).echelon_form().transpose().kernel().gens()[0]
	#	L2 = Matrix(Vectors).echelon_form().transpose().kernel().gens()[1]
	#	return [L1[i] + L2[i] for i in xrange(len(L1))]
	return list(Matrix(Vectors).echelon_form().transpose().kernel().gens()[0])

#-------------------------------------------------------------------------------
def GetNormalFromHNF(HNF):
	M = matrix(HNF)
	HNF = []
	for Row in M:
		HNF.append(Row)
	return FindNormal(HNF)

#-------------------------------------------------------------------------------
def DotProduct(U,V):
	DotProduct = 0
	for i in range(len(U)):
		DotProduct +=  U[i]*V[i]
	return DotProduct

#-------------------------------------------------------------------------------
def TransformPts(Pts, UCT):
	def ConvertMatrixToList(Matrix):
		Pts = []
		for Row in Matrix.transpose():
			Pts.append(list(Row))
		return Pts
	if len(Pts) < 1:
		return []
	TransformMatrix = UCT*(matrix(Pts).transpose())
	return ConvertMatrixToList(TransformMatrix)

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
def FindNewFacetPts(Pts, Edge, Normal, KnownFacetPts, NormalThroughFacet):
	FacetBarycenter = FindBarycenter(KnownFacetPts)
	MaxAngle = 'Test'
	NewFacetPts = []
	for Pt in Pts:
		if Pt not in KnownFacetPts:
			Vector = MakeVector(Edge[0], Pt)
			if n(DotProduct(Vector, Normal),1000) == 0:
				print "Denominator == 0!"
				raw_input()
			Num = DotProduct(Vector, NormalThroughFacet)
			Denom = DotProduct(Vector, Normal)
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
def FindNewFacetPtsFromEdge(Pts, Edge, Normal, KnownFacetPts):
	# Note that the input normal is an inner normal
	for i in xrange(len(Normal)):
		Normal[i] = Integer(Normal[i])				
	for Pt in Edge:
		for i in xrange(len(Pt)):
			Pt[i] = Integer(Pt[i])
	NormalThroughFacet = GetNormalFromHNF(GetHNFFromVectorAndPts(Normal, Edge))

	# I want this to be an outer facing normal
	for Pt in KnownFacetPts:
		if Pt not in Edge:
			AdditionalPtOnFacet = Pt
			break
	# For the sake of consistency, point the normal outwards. Decide later if that's what we want
	if NormalPointsTowardsPt(NormalThroughFacet, AdditionalPtOnFacet, Edge[0]):
		NormalThroughFacet = [-NormalThroughFacet[i] for i in xrange(len(NormalThroughFacet))]
	return FindNewFacetPts(Pts, Edge, Normal, KnownFacetPts, NormalThroughFacet)

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
		ValuesToPrepend = []
		for i in xrange(len(Row)):
			if Row[i] != "Diff":
				ValuesToPrepend.append(Row[i])
			if i != 0 and Row[i] != "Diff" and Row[i-1] == "Diff":
				print "Internal error: unexpected result from HNF in transformingpts"
				print Row
				print Pts
				raw_input()
		GoodPts = []
		for i in xrange(len(NewPts)):
			GoodPts.append([])
			for j in xrange(len(NewPts[i])):
				if Row[j] == "Diff":
					GoodPts[i].append(NewPts[i][j])
		ReverseUCT = matrix(UCT^-1, ZZ)
		PointMap = {}
		for Pt in GoodPts:
			PointMap[tuple(Pt)] = TransformPt(ValuesToPrepend + Pt, ReverseUCT)
	else:
		GoodPts = copy(Pts)
		PointMap = {}
		for Pt in GoodPts:
			PointMap[tuple(Pt)] = Pt
	return GoodPts, PointMap, Dim

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
			print "You found the points in an initial facet"
			return
	print "Pts you thought made a facet", PtsInFacet
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
		NewPts, Num, Denom = FindNewFacetPtsGood(Pts, FirstPts, Normal, FirstPts, NormalThroughFacet)
		Num, Denom = ScaleNumAndDenom(Num, Denom)
		FirstPts = NewPts + FirstPts
		PreviousNormal = copy(Normal)
		Normal = [-NormalThroughFacet[i]*Denom + Normal[i]*Num for i in xrange(len(Normal))]
		if not NormalPointsTowardsPt(Normal, Barycenter, FirstPts[0]):
			print "FLIP NORMAL"
			raw_input()
			Normal = [-Normal[i] for i in xrange(len(Normal))]
	CheckFacetAgainstSage(Pts, FirstPts)
	return FirstPts
