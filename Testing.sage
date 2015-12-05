load("RandomPretropisms.sage")

#-------------------------------------------------------------------------------
def DoTest(nvars, UseGiftWrap = False):
	PolyString = ""
	for i in xrange(nvars - 1):
		PolyString += "x_" + str(i) + ','
	PolyString += "x_" + str(nvars - 1)
	R = PolynomialRing(QQ, nvars, PolyString)
	HighestExp = 30
	NumberOfTerms = nvars + 1

	# I want a fast way to guarantee that I have full dimensional polytopes with the only nvars + 1 points.
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
	print "HI"

	global HullInfoMap
	HullInfoMap = {}
	HullTime = time()
	PtsList = []
	EdgeProd = 1

	for i in xrange(len(PolysAsPts)):
		print i, PolysAsPts[i]
		if UseGiftWrap:
			Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(PolysAsPts[i],True)
		else:
			Faces, IndexToPointMap, PointToIndexMap, Pts = FasterWrap(PolysAsPts[i])

		HullInfoMap[(i,"Faces")] = Faces
		HullInfoMap[(i,"IndexToPointMap")] = IndexToPointMap
		HullInfoMap[(i,"PointToIndexMap")] = PointToIndexMap
		HullInfoMap[(i,"Pts")] = Pts
		EdgeProd = EdgeProd*len(Faces[0])
		PtsList.append(Pts)

	print ""
	HullTime = time() - HullTime

	TimeList = []

	#DoGfan(Polys,R)
	TimeList.append(DoActualGfan(PolyString, GfanPolys))
	TempStr = ''
	#TimeList.append(Jeff(PolysAsPts, HullInfoMap,nvars - 1))
	TimeList.append(DoNewNewAlgorithm(PolysAsPts, HullInfoMap,nvars - 1)+HullTime)
	print "ConvexHullTime", HullTime
	#TimeList.append(DoNewAlgorithm(PolysAsPts, HullInfoMap))
	#TimeList.append(DoMinkowskiSum(Polys, PtsList))
	#TimeList.append(DoCayleyPolytope(PtsList))
	#TimeList.append(DoNaiveAlgorithm(HullInfoMap, nvars - 1))
	for i in xrange(6):
		print ""
	print TimeList
	return TimeList

#-------------------------------------------------------------------------------
def DoLotsOfTests(nvars, ntests):
	TempList = DoTest(nvars)
	for i in xrange(ntests-1):
		Trial = DoTest(nvars)
		TempList = [TempList[j] + Trial[j] for j in xrange(len(Trial))]
		print "WorkingSum of ", i + 2, "trials", TempList
	print "final result = ", [TempList[i]/ntests for i in xrange(len(TempList))]
	return [TempList[i]/ntests for i in xrange(len(TempList))]
