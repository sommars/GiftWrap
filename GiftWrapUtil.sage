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
	def RowIsZero(Row):
		for Int in Row:
			if Int != 0:
				return False
		return True
	Dim = 0
	for i in xrange(len(HNF[0])):
		if HNF[Dim][i] != 0:
			Dim += 1
			if Dim == len(HNF):
				break
	return Dim - 1

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
		return "AAASDFASFDASDFKASDF"
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
def NormVector(U):
	Norm = 0
	for i in range(len(U)):
		Norm += U[i]*U[i]
	return Norm**.5

#-------------------------------------------------------------------------------
def FindCosTheta(U,V):
#This feels wrong
	return float(DotProduct(U,V))/float((NormVector(U)*NormVector(V)))

#-------------------------------------------------------------------------------
def GetAMaximalPt(Pts, FacetPts, Normal):
	# We need some value to start with for theta, so pick the first valid option
	for i in xrange(len(Pts)):
		Pt = Pts[i]
		if (Pt not in FacetPts and not PtIsZero(Pt)):
			Theta = FindCosTheta(Pt, Normal)
			MaximalPt = Pt
			break

	# We don't need to go through all values, since we've already gone through some
	for j in xrange(i+1,len(Pts)):
		Pt = Pts[j]
		if (Pt not in FacetPts and not PtIsZero(Pt)):
			TestTheta = FindCosTheta(Pt, Normal)
			if TestTheta > Theta:
				MaximalPt = Pt
				Theta = TestTheta
	return MaximalPt

#-------------------------------------------------------------------------------
def NormalShouldBePositive(Pts, FacetPts):
	FacetVertValue = FacetPts[0][0]
	for i in xrange(1,len(FacetPts)):
		if FacetVertValue != FacetPts[i][0]:
			print "All of the facetpts aren't on the same hyperplane"
			raw_input()
			return "All of the facetpts aren't on the same hyperplane"
	# Sentinel value so that we know the value hasn't been set
	NormalShouldBePositive = -1
	for i in xrange(len(Pts)):
		if FacetVertValue > Pts[i][0]:
			if NormalShouldBePositive == -1:
				NormalShouldBePositive = True
			elif NormalShouldBePositive == False:
				print "All of the points aren't to the same side of the hyperplane"
				raw_input()
				return "All of the points aren't to the same side of the hyperplane"
		elif FacetVertValue < Pts[i][0]:
			if NormalShouldBePositive == -1:
				NormalShouldBePositive = False
			elif NormalShouldBePositive == True:	
				print "All of the points aren't to the same side of the hyperplane"
				raw_input()
				return "All of the points aren't to the same side of the hyperplane"
	return NormalShouldBePositive

#-------------------------------------------------------------------------------
def GetMaxAndMinPts(Pts, FacetPts, Normal):
	MaximalPts = []
	MinimalPts = []
	# We need some value to start with for the thetas, so pick the first valid option
	for i in xrange(len(Pts)):
		Pt = Pts[i]
		if (Pt not in FacetPts and not PtIsZero(Pt)):
			MaxTheta = FindCosTheta(Pt, Normal)
			MinTheta = MaxTheta
			MaximalPts.append(Pt)
			MinimalPts.append(Pt)
			break
	print Pt, i, MaxTheta
	# We don't need to go through all values, since we've already gone through some
	for j in xrange(i+1,len(Pts)):
		Pt = Pts[j]
		if (Pt not in FacetPts and not PtIsZero(Pt)):
			TestTheta = FindCosTheta(Pt, Normal)
			if TestTheta == MaxTheta:
				MaximalPts.append(Pt)
			if TestTheta == MinTheta:
				MinimalPts.append(Pt)
			if TestTheta > MaxTheta:
				MaximalPts = [Pt]
				MaxTheta = TestTheta
			if TestTheta < MinTheta:
				MinimalPts = [Pt]
				MinTheta = TestTheta
	return MaximalPts, MinimalPts

#-------------------------------------------------------------------------------
def PtIsZero(Pt):
	for Coord in Pt:
		if Coord != 0:
			return False
	return True

#-------------------------------------------------------------------------------
def TransformPts(Pts, UCT):
	def ConvertMatrixToList(Matrix):
		Pts = []
		for Row in Matrix.transpose():
			Pts.append(list(Row))
		return Pts
	TransformMatrix = UCT*(matrix(Pts).transpose())
	return ConvertMatrixToList(TransformMatrix)

#-------------------------------------------------------------------------------
def TransformPt(Pt, UCT):
	return TransformPts([Pt], UCT)[0]
