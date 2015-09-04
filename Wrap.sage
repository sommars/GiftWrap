load("GiftUtil.sage")
# The code in this unit is the code that will have to change for going between 3 and n dimensions.
def GiftWrap(Pts):
	SageHull = Polyhedron(vertices = Pts)
	Facets = []
	if not PtsAreValid(Pts):
		return "The input is not valid."
	Pts = RemoveDups(Pts)
	global Barycenter
	Barycenter = FindBarycenter(Pts)
	# Currently, I'm going to blow up if my pointset isn't full dimensional
	if CheckNormalFormDim(GetHNF(Pts)) != 2:
		print "NOT FULL DIMENSITONAL"
		print Pts
		return
	InitialFacet, PtsToRemove, EdgeList = FindInitialFacet(Pts)
	Facets.append(InitialFacet)
	Pts = RemovePts(Pts, PtsToRemove)
	# Going to spin through all of the edges in edgelist.
	# EdgeList should really be done via sets.
	Counter = 0
	while len(EdgeList) > 0:
		Edge = EdgeList[0]
		#Find Facet that has this Edge
		
		for Facet in Facets:
			Vertices = Facet.Vertices
			if (Edge[0] in Vertices) and (Edge[1] in Vertices):
				break
		
		Counter += 1
		if Counter == 200:
			return
			raw_input()
		FacetPts = Edge + FindNewFacetPtsFromEdgeTwo(Pts, Edge, Facet.InnerNormal, Facet.Vertices)

		#Call MakeFacet on the points I've found and do some bookkeeping
		Facet, PtsToRemove, Edges = MakeFacet(FacetPts, Pts)

		Pts = RemovePts(Pts, PtsToRemove)

		for Edge in Edges:
			# This is hacky. Really should just be careful and do some sorting
			ReversedEdge = [Edge[1], Edge[0]]
			if Edge in EdgeList:
				EdgeList.remove(Edge)
			elif ReversedEdge in EdgeList:
				EdgeList.remove(ReversedEdge)
			else:
				EdgeList.append(Edge)
		Facets.append(Facet)
	# If there are any points sitting in the interior of our polytope then they
	# have not yet been removed. We do so now. This is lazily programmed, and will
	# be fixed in the future.
	InteriorPts = list(Pts)
	for Facet in Facets:
		for Vertex in Facet.Vertices:
			if Vertex in InteriorPts:
				InteriorPts.remove(Vertex)

	Pts = RemovePts(Pts, InteriorPts)
	print "Finished! Count of faces = ", len(Facets)
	print ""
	print "....Computing Euler Characteristic...."
	print "Number of Vertices = ", len(Pts)
	EdgeCount = 0
	for Facet in Facets:
		#print Facet.Vertices
		EdgeCount += len(Facet.Vertices)
	EdgeCount = EdgeCount/2
	print "Number of Edges = ", EdgeCount
	print "Number of Facets", len(Facets)
	print len(Pts) - EdgeCount + len(Facets) == 2
	if len(Pts) - EdgeCount + len(Facets) != 2:
		print "Euler Characteristic failed!"
		if len(SageHull.faces(0)) != len(Pts):
			print "SageHullFaces0Length = ", len(SageHull.faces(0))
		if len(SageHull.faces(1)) != EdgeCount:
			print "SageHullFaces1Length = ", len(SageHull.faces(1))
		if len(SageHull.faces(2)) != len(Facets):
			print "SageHullFaces2Length = ", len(SageHull.faces(2))		
		raw_input()
	return Facets

#-------------------------------------------------------------------------------
def FindInitialFacet(Pts):
	print "Start initial facet"
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
	#print "hi!", FirstPts
	if len(FirstPts) == 1:
		#Start with one point
		FirstPts = FirstPts + FindNewFacetPtsFromSinglePt(Pts, FirstPts, [-1,0,0])
		#print "ASDF", FirstPts
		#Now we have either an edge or a facet. Need to check which
		if CheckNormalFormDim(GetHNF(FirstPts)) == 1:
			print "Final FirstPtsList: ", FirstPts
			return MakeFacet(FirstPts, Pts)

		KnownFacetPts = list(FirstPts)
		Counter = 1
		while True:
			NewPoint = [FirstPts[0][0],FirstPts[0][1]+Counter,FirstPts[0][2]]
			if not NewPoint in KnownFacetPts:
				KnownFacetPts.append(NewPoint)
				break
			Counter += 1
		#if we get here, we have an edge. need to pick up one or more points
		#First, we need to have an inner normal
		InnerNormal = GetNormalFromHNF(GetHNF(KnownFacetPts))
		"""
		print KnownFacetPts
		#CODE TO DELETE
		UCT = GetUCT(InnerNormal)
		NewPts = TransformPts(KnownFacetPts, UCT)
		VertHyperplaneX = NewPts[0][0]
		for i in xrange(1,len(NewPts)):
			if NewPts[i][0] != VertHyperplaneX:
				print NewPts
				raw_input()
		print UCT
		NewAllPts = TransformPts(Pts, UCT)
		VertHyperplaneIsGreatest = "maybe"
		for Pt in NewAllPts:
			if VertHyperplaneIsGreatest == True and Pt[0] > VertHyperplaneX:
				print VertHyperplaneX
				print NewAllPts
				raw_input()
			if VertHyperplaneIsGreatest == False and Pt[0] < VertHyperplaneX:
				print VertHyperplaneX
				print NewAllPts
				raw_input()
			if VertHyperplaneIsGreatest == "maybe":
				if Pt[0] < VertHyperplaneX:
					VertHyperplaneIsGreatest = True
				if Pt[0] > VertHyperplaneX:
					VertHyperplaneIsGreatest = False
		"""
		
		
		
		if not NormalPointsTowardsPt(InnerNormal, Barycenter, FirstPts[0]):
			InnerNormal = [-InnerNormal[i] for i in xrange(len(InnerNormal))]
		for i in xrange(len(InnerNormal)):
			InnerNormal[i] = Integer(InnerNormal[i])
		FirstPts = FirstPts + FindNewFacetPtsFromEdgeTwo(Pts, FirstPts, InnerNormal, KnownFacetPts)
		print "Final FirstPtsList: ", FirstPts
		return MakeFacet(FirstPts, Pts)		
	else:
		HNF = GetHNF(FirstPts)
		Dim = CheckNormalFormDim(HNF)
		if Dim == 1:
			print "Final FirstPtsList: ", FirstPts
			return MakeFacet(FirstPts, Pts)

#- However, if I have two points, this means that I have found an edge. I can do essentially the same thing I do in
#the general situation. I have the same normal to the hyperplane [-1,0,0]. I must use a specific normal to the pointset though, which I can #find by calling GetNormalFromHNF(GetHNFFromVectorAndPts([-1,0,0], Edge))
		elif Dim == 0:
			Normal = [-1,0,0]
			#Hacky
			KnownFacetPts = list(FirstPts)
			ExtraPt = [0 for i in xrange(len(FirstPts[0]))]
			ExtraPt[0] = FirstPts[0][0]
			for Pt in KnownFacetPts:
				for i in xrange(1,len(Pt)):
					ExtraPt[i] += abs(Pt[i])
			KnownFacetPts.append(ExtraPt)
			#EndHacky. Should refactor FindNewFacetPts routine so this isn't necessary
			FirstPts = FirstPts + FindNewFacetPtsFromEdgeTwo(Pts, FirstPts, [-1,0,0], KnownFacetPts)
			return MakeFacet(FirstPts, Pts)
		else:
			print "UNEXPECTED DIMENSION: ", Dim
			raw_input()
			return
	print "Final FirstPtsList: ", FirstPts
	return MakeFacet(FirstPts, Pts)

#-------------------------------------------------------------------------------
def MakeFacet(Pts, AllPts):
	# We have all of the points in their general position. We need to bring them
	# into vertical position so that the convexhull algorithm can work.
	UCT = GetUCTAndNormal(GetNormalFromHNF(GetHNF(Pts)))[0]
	NewPts = TransformPts(Pts, UCT)
	Hull, PointsToRemove = ConvexHull2d(NewPts)


	VertHyperplaneX = NewPts[0][0]
	for i in xrange(1,len(NewPts)):
		if NewPts[i][0] != VertHyperplaneX:
			print "Internal Error"
			print NewPts
			raw_input()

	MyAllPts = copy(AllPts)
	for Pt in Pts:
		MyAllPts.remove(Pt)
	NewAllPts = TransformPts(MyAllPts, UCT)
	VertHyperplaneIsGreatest = "maybe"
	for Pt in NewAllPts:
		if VertHyperplaneIsGreatest == True and Pt[0] > VertHyperplaneX:
			print VertHyperplaneX
			print NewAllPts
			raw_input()
		if VertHyperplaneIsGreatest == False and Pt[0] < VertHyperplaneX:
			print VertHyperplaneX
			print NewAllPts
			raw_input()
		if VertHyperplaneIsGreatest == "maybe":
			if Pt[0] < VertHyperplaneX:
				VertHyperplaneIsGreatest = True
			if Pt[0] > VertHyperplaneX:
				VertHyperplaneIsGreatest = False
		if Pt[0] == VertHyperplaneX:
			print "EQUALLL"
			print Pt
			print TransformPt(Pt, matrix(UCT^-1, ZZ))
			print TransformPts(Hull, matrix(UCT^-1, ZZ))
			print VertHyperplaneX
			print NewAllPts
			raw_input()
			
	Hull = TransformPts(Hull, matrix(UCT^-1, ZZ))
	
	
	#Hate doing this. Not sure why it's necessary, but this is a temporary solution
	for i in xrange(len(Hull)):
		for j in xrange(len(Hull[i])):
			Hull[i][j] = Integer(Hull[i][j])
	PointsToRemove = TransformPts(PointsToRemove, matrix(UCT^-1, ZZ))
	Facet = Face()

	Normal = GetNormalFromHNF(GetHNF(Hull))
	# Need to check if the normal is actually an innernormal
	if not NormalPointsTowardsPt(Normal, Barycenter, Hull[0]):
		Normal = [-Normal[i] for i in xrange(len(Normal))]
	Facet.InnerNormal = Normal
	Edges = [[Hull[len(Hull) - 1], Hull[0]]]
	for i in xrange(len(Hull) - 1):
		Edges.append([Hull[i],Hull[i+1]])
	#Hull.sort()
	Facet.Vertices = Hull
	Facet.Children = Hull
	Dim = CheckNormalFormDim(GetHNF(Hull))
	if Dim != 1:
		print "INCORRECT DIMENSION!!!"
		raw_input()
	Facet.Dimension = Dim
	Facet.Neighbors = []
	return Facet, PointsToRemove, Edges

	
for k in range(20):
	print ""
TestCube = [[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]
BigCube = []
for i in [0,1,2]:
	for j in [0,1,2]:
		for k in [0,1,2]:
			BigCube.append([i,j,k])
TetrahedronExtraPts = [[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1],[0,0,0],[1,0,0]]
Tetrahedron = [[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1]]
TippedOverHouse = [[0,0,0],[2,0,0],[0,2,0],[2,2,0],[0,0,2],[2,0,2],[0,2,2],[2,2,2],[3,1,1]]
Tests = [(TestCube,'Cube'), (Tetrahedron,'Tetrahedron'),(TippedOverHouse,'TippedOverHouse'),(TetrahedronExtraPts,'TetrahedronExtraPts'),(BigCube,'BigCube')]
#Tests = [(Tetrahedron,'Tetrahedron'),(TippedOverHouse,'TippedOverHouse'),(TetrahedronExtraPts,'TetrahedronExtraPts')]
#Tests = [(Tetrahedron,'Tetrahedron')]
#Tests = [(TetrahedronExtraPts,'TetrahedronExtraPts')]
#Tests = [(TippedOverHouse,'TippedOverHouse')]

def MakeRandomPointSet():
	Pts = []
	for i in xrange(randint(4,100)):
		Pts.append([])
		for j in xrange(3):
			Pts[i].append(Integer(randint(0,100)))
	return Pts
Tests = [(MakeRandomPointSet(), 'Random1')]


#Hang in main loop
#PASSED!
#ZZZZ = [[46, 48, 9], [99, 14, 6], [14, 91, 6], [77, 25, 81], [25, 72, 24], [56, 95, 46], [38, 2, 91], [56, 25, 79], [13, 96, 50], [43, 45, 22], [30, 89, 72]]

#PASSED!
#ZZZZ = [[98, 15, 16], [11, 10, 99], [73, 84, 41], [93, 91, 69], [67, 85, 52], [25, 24, 49], [62, 77, 65], [56, 75, 33], [40, 33, 55], [74, 53, 53]]

#PASSED!
ZZZZ = [[94, 28, 20], [92, 33, 93], [3, 81, 1], [73, 75, 15], [93, 90, 77]]

#Euler Characteristic failed
#PASSED!
#ZZZZ = [[7, 99, 74], [10, 31, 38], [29, 7, 23], [61, 7, 72], [29, 91, 1], [19, 8, 27], [66, 9, 27], [51, 50, 2], [32, 96, 6], [17, 59, 0]]

#Hangs indefinitely in the main loop
#Works with increased tolerance!
#ZZZZ = [[72, 19, 88], [74, 40, 53], [45, 44, 58], [43, 92, 13], [10, 6, 16], [50, 45, 95], [33, 6, 59], [84, 74, 58], [88, 91, 20], [1, 100, 69], [35, 62, 44], [27, 57, 25], [76, 9, 78], [15, 50, 38], [31, 91, 19], [29, 65, 33], [77, 73, 88], [98, 3, 48], [30, 16, 1], [50, 93, 65], [41, 62, 33], [35, 77, 77], [57, 4, 22], [61, 83, 99], [81, 53, 35], [73, 37, 74], [99, 24, 92], [8, 64, 48], [73, 36, 5], [59, 60, 99], [94, 51, 12], [74, 50, 72], [0, 65, 20], [57, 92, 3], [68, 71, 77], [86, 5, 17], [33, 58, 98], [58, 54, 51], [49, 46, 47], [26, 54, 82], [14, 6, 25], [3, 66, 27], [57, 81, 15], [73, 0, 61], [19, 21, 54], [97, 30, 50], [8, 27, 49], [50, 29, 90], [3, 92, 77], [55, 64, 54], [70, 85, 17], [56, 40, 43], [18, 31, 67], [33, 21, 99], [75, 82, 50], [75, 33, 50], [34, 8, 19], [97, 75, 13], [71, 60, 95], [26, 26, 2], [46, 1, 83], [26, 75, 91], [72, 75, 55], [55, 98, 23], [74, 4, 12], [30, 18, 28], [82, 72, 48], [55, 90, 36], [12, 85, 41], [67, 55, 21], [49, 97, 99], [90, 73, 9], [95, 84, 18], [95, 7, 81], [29, 86, 35], [62, 48, 76], [94, 7, 98]]


#Hangs indefinitely in the main loop
#Only broken because findinitialfacet is broken!
#ZZZZ = [[81, 47, 67], [95, 63, 1], [85, 44, 84], [51, 62, 58], [73, 45, 65], [84, 38, 42], [48, 72, 91], [45, 89, 64], [9, 70, 78], [0, 59, 89], [39, 11, 70], [21, 69, 99], [28, 10, 18], [59, 40, 39], [20, 45, 19], [21, 73, 15], [55, 77, 93], [76, 54, 25], [96, 1, 68], [70, 77, 15], [42, 59, 50], [89, 33, 79], [11, 80, 74], [65, 37, 100], [5, 100, 73], [52, 54, 0], [0, 78, 68], [50, 67, 18], [36, 78, 68], [99, 37, 71], [65, 13, 31], [48, 51, 45], [91, 48, 43], [96, 13, 22], [30, 18, 38], [86, 29, 30], [23, 88, 94], [90, 87, 48], [64, 34, 44], [80, 93, 14], [4, 67, 82], [37, 55, 52], [21, 44, 36], [91, 95, 79]]


#Float division by zero
#INITIALFACETFAILED!
# These two fail because i'm "randomly" picking the first normals. This is fixable
#ZZZZ = [[97, 74, 82], [35, 27, 85], [77, 94, 6], [92, 96, 46]]
#ZZZZ = [[82, 64, 49], [92, 97, 44], [82, 37, 76], [72, 21, 75]]

#Denom = 0
#Fixed!
#ZZZZ = [[90, 37, 29], [73, 45, 29], [59, 3, 61], [8, 52, 37], [8, 63, 62], [87, 51, 8], [1, 65, 8], [26, 44, 12], [16, 77, 74], [72, 50, 39], [55, 71, 57], [21, 25, 4], [4, 24, 22], [11, 100, 21], [22, 25, 88], [85, 23, 100], [20, 61, 16], [32, 54, 5], [38, 41, 93], [7, 92, 59], [89, 100, 63], [62, 85, 59], [28, 46, 54], [45, 64, 36], [70, 88, 12], [23, 46, 99], [61, 16, 96], [88, 59, 19], [62, 85, 92], [71, 36, 35], [32, 85, 13], [58, 89, 4], [64, 39, 36], [19, 35, 55], [21, 49, 73], [8, 87, 61], [38, 54, 1], [35, 86, 18], [53, 84, 69], [56, 73, 66], [86, 20, 90], [84, 96, 44], [13, 7, 90], [58, 70, 98], [68, 24, 41], [82, 24, 58], [18, 88, 7], [31, 56, 19], [20, 69, 17], [23, 0, 5], [43, 53, 5], [17, 65, 35], [42, 3, 8], [55, 77, 9], [89, 34, 34], [52, 28, 55], [100, 55, 72], [65, 74, 72], [42, 47, 35], [60, 51, 15], [54, 4, 45], [19, 4, 63], [94, 59, 41], [60, 35, 94], [18, 7, 83], [35, 0, 77], [8, 75, 2], [76, 33, 34]]
				

#Euler Failed...badly
#Fixed!
#ZZZZ = [[18, 14, 3], [100, 25, 29], [94, 81, 61], [37, 51, 84], [8, 74, 18], [31, 84, 67], [94, 96, 70], [72, 91, 0], [81, 39, 16], [48, 73, 85], [80, 19, 58], [93, 4, 30], [94, 90, 47], [32, 64, 89], [39, 88, 37], [78, 21, 42], [60, 84, 73], [71, 25, 37], [14, 68, 78], [51, 11, 75], [1, 6, 16], [66, 30, 73], [35, 82, 65], [49, 25, 65], [43, 22, 30], [27, 27, 71], [58, 8, 75], [17, 64, 65], [65, 11, 79], [2, 32, 61], [15, 5, 22], [53, 66, 86], [39, 81, 72], [10, 93, 58], [63, 24, 72], [98, 59, 53], [32, 36, 6], [6, 84, 40], [73, 68, 37], [51, 8, 53], [97, 29, 53], [51, 74, 60], [94, 33, 18], [58, 17, 24], [5, 10, 13], [49, 89, 13], [97, 66, 70], [99, 91, 46], [61, 57, 5], [41, 4, 77], [24, 2, 90], [51, 35, 73], [23, 62, 2], [50, 11, 29], [90, 75, 20], [74, 12, 13], [47, 75, 65], [51, 17, 84], [88, 81, 72], [66, 50, 30], [27, 29, 54], [38, 61, 15], [80, 38, 98], [98, 81, 33], [99, 31, 24], [57, 41, 42], [20, 79, 4], [32, 75, 84], [74, 72, 85], [19, 33, 64], [48, 94, 70], [15, 82, 68], [82, 68, 18], [18, 19, 55], [82, 5, 97], [95, 80, 25], [67, 51, 52], [97, 40, 44], [6, 69, 19]]

#Can't find initial facet
#Fixed! 
#ZZZZ = [[38, 11, 29], [23, 74, 62], [24, 7, 62], [21, 3, 29], [87, 74, 17], [1, 86, 18], [72, 60, 33], [19, 39, 65], [89, 81, 6], [50, 84, 16], [31, 22, 20], [38, 18, 54], [6, 96, 45], [74, 67, 80], [5, 95, 63], [1, 93, 29], [26, 87, 1], [48, 73, 41], [17, 97, 37], [69, 33, 1], [74, 72, 28], [65, 10, 30], [14, 28, 21], [71, 53, 99], [77, 58, 26], [85, 1, 63], [80, 63, 43], [40, 5, 58], [64, 29, 62], [37, 7, 2], [36, 25, 5], [6, 93, 55], [1, 6, 29], [35, 40, 13], [85, 71, 42], [93, 33, 16], [42, 27, 59], [71, 88, 15], [99, 50, 39], [46, 26, 1], [40, 64, 98], [78, 35, 28], [10, 91, 64], [84, 98, 68], [73, 22, 74], [83, 41, 42], [33, 24, 53], [63, 77, 14], [79, 51, 53], [24, 81, 87], [66, 27, 14], [82, 18, 47], [12, 12, 74], [98, 62, 71], [49, 57, 54], [42, 87, 19], [20, 19, 55], [89, 24, 73], [63, 14, 32], [6, 22, 40], [57, 73, 40], [24, 42, 65], [22, 61, 63], [89, 97, 40], [81, 70, 24], [56, 50, 84], [68, 80, 38], [39, 35, 97]]

#Euler char failed
#FIXED!
#ZZZZ = [[22, 88, 49], [24, 71, 64], [13, 78, 93], [56, 21, 98], [94, 16, 85], [35, 1, 39], [45, 82, 4], [87, 11, 78], [68, 31, 0], [20, 4, 31], [5, 81, 0], [78, 29, 48], [52, 65, 3], [83, 4, 13], [60, 71, 67], [85, 87, 40], [28, 100, 20], [34, 23, 10], [72, 20, 4], [11, 4, 98], [31, 36, 20], [43, 82, 79], [2, 47, 88], [12, 95, 8], [5, 39, 76], [1, 73, 83], [66, 87, 50], [25, 94, 12], [38, 75, 94], [75, 37, 48], [93, 83, 11], [51, 31, 12], [13, 43, 41], [65, 62, 14]]

#ZZZZ = [[1, 11, 45], [50, 87, 17], [61, 48, 50], [51, 64, 9], [3, 92, 36], [4, 78, 66], [90, 18, 64], [1, 22, 41], [14, 68, 81], [14, 14, 61], [92, 94, 63], [42, 80, 71], [55, 33, 93], [57, 80, 4], [36, 40, 32], [89, 40, 77], [15, 68, 60], [94, 98, 69]]

#ZZZZ = [[100, 31, 67], [68, 45, 7], [18, 44, 80], [97, 84, 74], [18, 46, 77], [5, 95, 57], [99, 42, 1], [84, 35, 75], [62, 41, 33], [62, 95, 47], [30, 93, 7], [28, 28, 85], [38, 23, 61], [44, 28, 28], [47, 54, 39], [34, 16, 10], [35, 80, 85], [49, 21, 44], [66, 21, 5], [41, 87, 56], [94, 90, 94], [27, 24, 85], [81, 26, 18], [30, 31, 53], [75, 68, 46], [95, 14, 43], [40, 34, 46], [50, 1, 48], [84, 38, 93], [9, 85, 96], [21, 85, 37], [94, 60, 8], [51, 91, 82], [82, 96, 46], [23, 23, 41], [11, 65, 100], [2, 1, 10], [1, 66, 44], [62, 83, 65], [52, 9, 20], [70, 41, 84], [14, 7, 15], [83, 72, 42], [63, 57, 40], [17, 100, 76], [98, 67, 68], [52, 62, 44], [22, 21, 92], [88, 96, 85], [58, 98, 34], [24, 94, 38], [5, 13, 67], [94, 75, 97]]

#ZZZZ = [[76, 66, 23], [35, 44, 47], [25, 3, 36], [96, 35, 88], [74, 22, 62], [46, 13, 71], [50, 80, 33], [48, 12, 79], [21, 75, 58], [34, 40, 20], [41, 13, 76], [41, 58, 3], [31, 21, 22], [6, 12, 74], [79, 25, 29], [62, 27, 15], [66, 88, 25], [27, 0, 9], [0, 26, 60], [28, 93, 99], [58, 36, 2], [7, 87, 95], [66, 100, 7], [62, 93, 13], [13, 22, 9], [48, 91, 20], [87, 94, 23], [90, 3, 29], [62, 13, 77], [92, 85, 91], [18, 38, 46], [84, 34, 69], [25, 82, 22], [22, 39, 97], [32, 84, 19], [61, 59, 40], [1, 20, 30], [31, 8, 42], [84, 58, 41], [97, 68, 94], [69, 86, 53], [78, 20, 54], [80, 72, 16], [26, 2, 5], [0, 30, 32], [97, 45, 71], [93, 34, 49], [11, 99, 23], [93, 98, 23], [16, 2, 26], [49, 81, 61], [20, 58, 20], [86, 6, 89], [19, 62, 99], [23, 91, 92], [54, 11, 10], [43, 67, 23], [99, 1, 66], [14, 15, 19], [6, 35, 69], [83, 32, 21], [81, 38, 44], [62, 40, 15], [51, 5, 71], [83, 62, 81], [28, 47, 73], [86, 43, 88], [5, 54, 79], [8, 78, 77], [32, 77, 78], [67, 13, 17], [76, 81, 75], [15, 80, 91], [37, 97, 18], [47, 66, 53], [61, 99, 38], [49, 62, 56], [16, 32, 84], [40, 45, 32], [82, 77, 81], [100, 88, 92], [59, 53, 84], [10, 7, 66], [19, 46, 56], [9, 98, 70], [4, 71, 24], [2, 64, 16], [95, 69, 78], [79, 50, 58], [43, 16, 48], [94, 45, 42], [53, 73, 84]]


#Denom = 0
#For some reason, this facet: [[89, 100, 94], [97, 58, 0], [92, 93, 85]] lies on the same hyperplane as this point: [78, 94, 32]
# not yet sure why this is happening....
ZZZZ = [[74, 58, 61], [57, 100, 36], [46, 68, 77], [2, 23, 2], [62, 68, 57], [69, 23, 69], [37, 3, 34], [3, 74, 97], [37, 62, 69], [65, 65, 84], [47, 48, 60], [20, 78, 58], [73, 39, 63], [95, 10, 50], [24, 63, 66], [20, 49, 75], [10, 95, 86], [25, 79, 6], [69, 17, 53], [83, 55, 72], [51, 72, 83], [100, 50, 40], [30, 29, 31], [24, 49, 39], [1, 80, 98], [4, 42, 80], [29, 24, 20], [98, 26, 69], [98, 22, 65], [12, 23, 53], [89, 100, 94], [39, 73, 30], [50, 1, 3], [30, 82, 12], [64, 50, 5], [73, 89, 5], [12, 76, 40], [87, 53, 97], [5, 89, 21], [71, 19, 53], [33, 51, 36], [16, 3, 31], [91, 6, 2], [50, 48, 92], [52, 82, 39], [53, 84, 60], [50, 39, 17], [63, 52, 34], [76, 34, 40], [39, 11, 41], [94, 93, 98], [33, 70, 92], [8, 53, 12], [88, 27, 19], [97, 58, 0], [75, 21, 33], [49, 87, 35], [74, 81, 56], [73, 66, 54], [29, 97, 60], [77, 75, 4], [86, 33, 66], [75, 99, 93], [0, 65, 34], [60, 66, 26], [63, 69, 45], [16, 52, 53], [74, 96, 38], [78, 94, 32], [92, 55, 57], [27, 42, 60], [1, 56, 83], [12, 21, 73], [80, 56, 50], [29, 77, 80], [56, 3, 29], [6, 37, 14], [55, 17, 83], [41, 48, 67], [70, 95, 45], [72, 91, 60], [92, 93, 85], [70, 81, 75], [40, 51, 27]]

#same..
ZZZZ = [[50, 80, 92], [28, 0, 44], [47, 87, 12], [3, 41, 11], [53, 22, 9], [85, 100, 78], [15, 24, 18], [31, 29, 52], [11, 23, 19], [86, 66, 89], [80, 50, 56], [97, 18, 46], [16, 9, 36], [44, 20, 24], [11, 51, 68], [7, 9, 69], [15, 97, 45], [55, 78, 100], [13, 64, 27], [6, 64, 78], [30, 41, 79], [53, 30, 99], [15, 11, 7], [99, 29, 25], [0, 0, 58], [10, 98, 83], [30, 54, 75], [39, 14, 100], [85, 38, 49], [2, 70, 62], [29, 48, 48], [58, 86, 32], [10, 25, 67], [34, 8, 47], [79, 69, 73], [21, 96, 77], [31, 76, 34], [87, 22, 55], [48, 17, 41], [23, 69, 9], [95, 10, 58], [50, 13, 42], [62, 75, 1], [74, 59, 79], [76, 58, 26], [40, 10, 100], [82, 12, 52], [15, 51, 80], [31, 48, 90], [79, 2, 90], [50, 27, 52], [49, 21, 66], [39, 32, 52], [17, 60, 66], [98, 77, 13], [26, 99, 63], [34, 43, 18], [100, 15, 76], [8, 43, 67], [21, 26, 39], [48, 54, 98], [19, 5, 94], [54, 33, 87], [51, 36, 93], [27, 73, 76], [92, 12, 4], [88, 68, 8], [54, 6, 80], [72, 92, 79], [12, 64, 89], [27, 33, 9], [3, 66, 17], [14, 55, 80], [52, 73, 93], [15, 32, 29], [88, 79, 51], [18, 18, 69], [64, 23, 64], [28, 30, 100], [86, 7, 75], [94, 21, 86], [12, 31, 59], [96, 99, 4], [30, 61, 80], [22, 89, 61], [54, 73, 61], [94, 4, 28], [84, 35, 30], [63, 37, 74], [16, 18, 92], [81, 73, 73], [20, 61, 1], [19, 32, 47], [19, 67, 88], [58, 41, 99], [65, 50, 80], [28, 100, 42], [55, 30, 57], [61, 70, 100], [83, 43, 62]]

ZZZZ = [[17, 0, 49], [69, 39, 27], [3, 64, 8], [94, 82, 96], [81, 100, 50], [60, 14, 29], [0, 95, 10], [40, 36, 76], [33, 17, 16], [25, 21, 71], [97, 27, 5], [5, 29, 54], [2, 98, 30], [30, 57, 51], [29, 43, 82], [60, 20, 100], [58, 38, 25], [31, 66, 66], [21, 60, 27], [1, 54, 31], [70, 97, 40], [34, 25, 87], [95, 14, 36], [87, 73, 1], [25, 60, 4], [29, 92, 99], [91, 26, 95], [88, 87, 17], [50, 60, 17], [41, 100, 69], [58, 45, 23], [92, 18, 99], [83, 1, 95], [4, 36, 71], [84, 69, 98], [9, 50, 33], [99, 50, 44], [18, 28, 6], [90, 22, 89], [84, 49, 22], [81, 15, 96], [13, 0, 16], [38, 79, 87], [61, 31, 56], [62, 45, 35], [56, 5, 35], [9, 59, 51], [100, 76, 7], [85, 0, 3], [26, 85, 53], [77, 51, 5], [63, 100, 23], [37, 4, 86], [84, 91, 62], [26, 1, 96], [95, 0, 28], [27, 24, 97], [25, 83, 4], [16, 72, 24], [80, 88, 65], [92, 76, 85], [99, 68, 15], [66, 22, 31]]

ZZZZ = [[50, 89, 83], [34, 61, 76], [92, 15, 80], [43, 91, 70], [89, 43, 28], [62, 50, 72], [80, 5, 66], [44, 24, 80], [27, 18, 30], [64, 38, 89], [42, 39, 16], [73, 97, 95], [91, 8, 90], [30, 67, 24], [0, 72, 79], [24, 13, 15], [66, 29, 36], [4, 47, 93], [46, 2, 68], [24, 26, 13], [2, 78, 69], [73, 19, 13], [84, 98, 29], [64, 32, 4], [80, 23, 42], [95, 87, 32], [63, 15, 37], [94, 19, 28], [98, 34, 50], [89, 79, 88], [75, 13, 78], [95, 23, 50], [82, 50, 63], [49, 70, 32], [2, 24, 69], [62, 85, 14], [82, 88, 7], [54, 79, 88], [72, 31, 77], [15, 41, 46], [91, 28, 94], [62, 96, 63], [91, 31, 94], [0, 87, 80], [24, 77, 70], [100, 24, 22], [68, 48, 31], [13, 100, 38], [38, 88, 11], [92, 0, 24], [47, 36, 7], [98, 60, 87], [91, 98, 43], [7, 44, 51], [42, 62, 58], [40, 89, 94], [59, 30, 47], [30, 69, 94], [56, 35, 32], [89, 60, 96], [36, 99, 58]]
	
		
Tests = []
for i in xrange(500):
	Tests.append((MakeRandomPointSet(), 'Random1'))
	
#Tests = [(ZZZZ, 'ZZZZ')]

#Tests = [(TestCube,'Cube'), (Tetrahedron,'Tetrahedron'),(TippedOverHouse,'TippedOverHouse'),(TetrahedronExtraPts,'TetrahedronExtraPts'),(BigCube,'BigCube')]
#Tests = [(TippedOverHouse,'TippedOverHouse')]

#Tests = [(TestCube,'Cube')]

for Test in Tests:
	Pts = Test[0]
	print Test[1]
	print Pts
#	print FindInitialFacet(Pts)
	GiftWrap(Pts)
	print "Done with " + Test[1]
	print ""
	print ""
