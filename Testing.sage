load("RandomPretropisms.sage")

#-------------------------------------------------------------------------------
def DoReducedCyclicTest(nvars, UseGiftWrap = False):
	"""
	Function that allows speed comparisons of many different methods to compute
	pretropisms for the cyclic-n problem.
	"""
	
	Polys, R = ReducedCyclicNRoots(nvars)
	#Manually reordering the polynomials by size of lineality space
	if nvars == 5:
		Polys = Polys
	if nvars == 6:
		Polys = [Polys[0],Polys[4],Polys[1],Polys[3],Polys[2]]
	if nvars == 7:
		Polys = Polys
	if nvars == 8:
		Polys = [Polys[0],Polys[2],Polys[4],Polys[6],Polys[1],Polys[5],Polys[3]]
	if nvars == 9:
		Polys = [Polys[0],Polys[1],Polys[3],Polys[4],Polys[6],Polys[7],Polys[2],Polys[5]]
	if nvars == 10:
		Polys = [Polys[0],Polys[2],Polys[6],Polys[8],Polys[1],Polys[3],Polys[5],Polys[7],Polys[4]]
	if nvars == 11:
		Polys = Polys	
	if nvars == 12:
		Polys = [Polys[0],Polys[4],Polys[6],Polys[10],Polys[1],Polys[9],Polys[2],Polys[8],Polys[3],Polys[7],Polys[5]]
	
	PolysAsPts = [[[Integer(j) for j in i] for i in Poly.exponents()] for Poly in Polys]

	global HullInfoMap
	HullInfoMap = {}
	global NumberOfPolytopes
	NumberOfPolytopes = len(PolysAsPts)
	HullTime = time()
	for i in xrange(len(PolysAsPts)):
		Faces, IndexToPointMap, PointToIndexMap, Pts = WrapWithLineality(PolysAsPts[i], R, [Polys[i]])
		HullInfoMap[(i,"Faces")] = Faces
		HullInfoMap[(i,"IndexToPointMap")] = IndexToPointMap
		HullInfoMap[(i,"PointToIndexMap")] = PointToIndexMap
		HullInfoMap[(i,"Pts")] = Pts
	HullTime = time() - HullTime
	print "ConvexHullTime", HullTime
	TimeList = [HullTime]
	TimeList.append(DoGfanFromSage(Polys, R))
	TimeList.append(DoNewAlgorithm(HullInfoMap))
	#TimeList.append(DoNewAlgorithmTwo(HullInfoMap,[i for i in xrange(len(PolysAsPts))]))
	print TimeList
	return TimeList

#-------------------------------------------------------------------------------
def DoGenericTest(nvars, UseGiftWrap = False):
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
	print PolysAsPts
	#Instead of having a lot of lists sitting around, I have the face structures,
	#maps, and a list of the pts all stored in a map. Each polytope has its values
	#stored here
	global HullInfoMap
	global NumberOfPolytopes
	NumberOfPolytopes = len(PolysAsPts)
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
	TimeList.append(DoNewAlgorithm(HullInfoMap))
	#TimeList.append(DoMinkowskiSum(Polys, PtsList))
	#TimeList.append(DoCayleyPolytope(PtsList))
	#TimeList.append(DoGfan(PolyString, GfanPolys))
	TimeList.append(DoGfanFromSage(Polys, R))
	#TimeList.append(DoConeIntersectionAlgorithm(HullInfoMap, nvars - 1))

	print TimeList
	return TimeList

#-------------------------------------------------------------------------------
def DoLotsOfTests(nvars, ntests):
	"""
	This is a wrapper so that DoGenericTest can easily be called many times.
	"""
	if ntests < 1 or nvars < 3:
		print "Invalid input."
		return
	TempList = DoGenericTest(nvars)
	for i in xrange(1,ntests):
		Trial = DoGenericTest(nvars)
		TempList = [TempList[j] + Trial[j] for j in xrange(len(Trial))]
		print "WorkingSum of ", i + 1, "trials", TempList
	print "final result = ", [TempList[i]/ntests for i in xrange(len(TempList))]
	return [TempList[i]/ntests for i in xrange(len(TempList))]
