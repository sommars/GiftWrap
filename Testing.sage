load("RandomPretropisms.sage")

#-------------------------------------------------------------------------------
def DoTest(nvars, UseGiftWrap = False):
	"""
	Function that allows speed comparisons of many different methods to compute
	pretropisms.
	"""
	
	PolyString = ""
	for i in xrange(nvars - 1):
		PolyString += "x_" + str(i) + ','
	PolyString += "x_" + str(nvars - 1)
	R = PolynomialRing(QQ, nvars, PolyString)
	HighestExp = 30
	NumberOfTerms = nvars + 1

	#I want a fast way to guarantee that I have full dimensional simplicial polytopes
	Polys = []
	GfanPolys = []
	while len(Polys) < nvars - 1:
		TestPoly = R.random_element(HighestExp,NumberOfTerms)
		if Polyhedron([[Integer(j) for j in i] for i in TestPoly.exponents()]).dim() == nvars:
			Polys.append(TestPoly)
			NewGfpoly = 0
			for Monomial in TestPoly.monomials():
				NewGfpoly += Monomial
			GfanPolys.append(NewGfpoly)

	PolysAsPts = [[[Integer(j) for j in i] for i in Poly.exponents()] for Poly in Polys]

	#Instead of having a lot of lists sitting around, I have the face structures,
	#maps, and a list of the pts all stored in a map. Each polytope has its values
	#stored here
	global HullInfoMap
	HullInfoMap = {}
	HullTime = time()
	for i in xrange(len(PolysAsPts)):
		if UseGiftWrap:
			Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(PolysAsPts[i],True)
		else:
			Faces, IndexToPointMap, PointToIndexMap, Pts = FasterWrap(PolysAsPts[i])
		HullInfoMap[(i,"Faces")] = Faces
		HullInfoMap[(i,"IndexToPointMap")] = IndexToPointMap
		HullInfoMap[(i,"PointToIndexMap")] = PointToIndexMap
		HullInfoMap[(i,"Pts")] = Pts

	HullTime = time() - HullTime
	print "ConvexHullTime", HullTime
	TimeList = [HullTime]
	TimeList.append(DoNewAlgorithm(PolysAsPts, HullInfoMap,nvars - 1))
	#TimeList.append(DoMinkowskiSum(Polys, PtsList))
	#TimeList.append(DoCayleyPolytope(PtsList))
	#TimeList.append(DoGfan(PolyString, GfanPolys))
	#TimeList.append(DoGfanFromSage(Polys, R))
	#TimeList.append(DoConeIntersectionAlgorithm(HullInfoMap, nvars - 1))

	print TimeList
	return TimeList

#-------------------------------------------------------------------------------
def DoLotsOfTests(nvars, ntests):
	"""
	This is a wrapper so that DoTest can easily be called many times.
	"""
	if ntests < 1 or nvars < 3:
		print "Invalid input."
		return
	TempList = DoTest(nvars)
	for i in xrange(1,ntests):
		Trial = DoTest(nvars)
		TempList = [TempList[j] + Trial[j] for j in xrange(len(Trial))]
		print "WorkingSum of ", i + 1, "trials", TempList
	print "final result = ", [TempList[i]/ntests for i in xrange(len(TempList))]
	return [TempList[i]/ntests for i in xrange(len(TempList))]
