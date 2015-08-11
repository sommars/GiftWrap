load("GiftUtil.sage")
# The code in this unit is the code that will have to change for going between 3 and n dimensions.
def GiftWrap(Pts):
	Facets = []
	if not PtsAreValid(Pts):
		return "The input is not valid."
	Pts = RemoveDups(Pts)
	InitialFacet, PtsToRemove, EdgeList = FindInitialFacet(Pts)
	Facets.append(InitialFacet)
	Pts = RemovePts(Pts, PtsToRemove)
	# Going to spin through all of the edges in edgelist.
	# EdgeList should really be done via sets.
	while len(EdgeList) > 0:
		Edge = EdgeList[0]
		#Find Facet that has this Edge
		for Facet in Facets:
			Vertices = Facet.Vertices
			if (Edge[0] in Vertices) and (Edge[1] in Vertices):
				break
		#Do a UCT on the normal to this edge
		UCT, Normal = GetUCTAndNormal(GetNormalFromHNF(GetHNF(Edge)))
		NewEdge = TransformPts(Edge,UCT)
		#Do the same UCT on the pointset
		NewPts = TransformPts(Pts, UCT)
		#Get the additional points on this facet (Facet points will be the edge plus whatever we find in the GetMaximalPts)
		MaximalPts, MinimalPts = GetMaxAndMinPts(NewPts, NewEdge, Normal)
		# We have both maximal and minimal points. One direction is the way we want
		# to rotate, the other will give us the facet that had that edge previously.
		# We check which points we want here.
		ExistingFacetVertices = TransformPts(Vertices, UCT)
		UseMaximalPts = False
		for Pt in MaximalPts:
			if Pt in ExistingFacetVertices:
				UseMaximalPts = True
		if UseMaximalPts == True:
			AdditionalPts = MaximalPts
		else:
			AdditionalPts = MinimalPts

		FacetPts = Edge + TransformPts(AdditionalPts, matrix(UCT^-1,ZZ))
		# I may have found one or more additional pts on the hyperplane. I don't know
		# that I've found all of them. Do a UCT on the FacetPts, see what else is in
		# the hyperplane
		UCT = GetUCT(FacetPts)
		NewPts = TransformPts(Pts, UCT)
		NewFacetPts = TransformPts(FacetPts, UCT)
		for Pt in NewPts:
			if (Pt[0] == 1) and (Pt not in NewFacetPts):
				NewFacetPts.append(Pt)
		#Call MakeFacet on the points I've found and do some bookkeeping
		Facet, PtsToRemove, Edges = MakeFacet(FacetPts)
		Pts = RemovePts(Pts, PtsToRemove)
		for Edge in Edges:
			if Edge in EdgeList:
				EdgeList.remove(Edge)
			else:
				EdgeList.append(Edge)
		Facets.append(Facet)
	for Facet in Facets:
		Facet.PrintProps()
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
	print FirstPts
	Counter = 0
	while True:
		Counter += 1
		HNF = GetHNF(FirstPts)
		if CheckNormalFormDim(HNF) == 2:
			break
		Normal = GetNormalFromHNF(HNF)
		# If we're going through this the first time and we have a zero dimensional
		# set of points found so far, we want the normal to be the vector which we
		# minimized, i.e. [1,0,...,0], rather than what we get by taking the Hermite
		# Normal Form normal.
		if (Counter == 1) or len(FirstPts) == 1:	
			Normal[0] = 1
			for i in xrange(1,len(Normal)):
				Normal[i] = 0
		UCT, NewNormal = GetUCTAndNormal(Normal)
		print ""
		print "ROUND: ", Counter
		print "Normal:" , Normal
		print HNF
		print "UCTOFNORMAL:", TransformPt(NewNormal[0], UCT)
		print "NewNormal:" , NewNormal
		print "UCT:"
		print UCT
		print "FirstPts", FirstPts
		print "TransformedFirstPts", TransformPts(FirstPts, UCT)
		print "ALLPTSTRANSFORMED", TransformPts(Pts, UCT)
		if len(NewNormal) > 1:
			return "This is a disaster!!!!"
		NewNormal = NewNormal[0]
		NewPts = TransformPts(Pts, UCT)
		# Find the new points and add them to FirstPts. We add the points in the
		# original space since FirstPts is in the original space. Also, we transform
		# the normal to the new space, because otherwise its value is meaningless.
		FirstPts = FirstPts + TransformPts(GetMaximalPts(NewPts, FirstPts, NewNormal), matrix(UCT^-1,ZZ))
	print FirstPts
	return MakeFacet(FirstPts)
#Question for Jan: why is it that the normal doesn't also transform?

#-------------------------------------------------------------------------------
def MakeFacet(Pts):
	Hull, PointsToRemove = ConvexHull2d(Pts)
	Facet = Face()
#	Facet.Normal = 
	Facet.Vertices = Hull
	Facet.Children = Hull
	Facet.Dimension = 2
	Facet.Neighbors = []
	Edges = [[Hull[len(Hull) - 1], Hull[0]]]
	for i in xrange(len(Hull) - 1):
		Edges.append([Hull[i],Hull[i+1]])
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
Tests = [(Tetrahedron,'Tetrahedron'),(TippedOverHouse,'TippedOverHouse'),(TetrahedronExtraPts,'TetrahedronExtraPts')]


Tests = [(TippedOverHouse,'TippedOverHouse')]
for Test in Tests:
	Pts = Test[0]
	print Test[1]
	print Pts
	#GiftWrap(Pts)
	FindInitialFacet(Pts)
	print "Done with " + Test[1]
	"""
	print "Number of Vertices"
	#Note that this needs to be the actual vertices not the start set of pts
	print len(Pts)
	print "Number of Edges"
	print "Number of Facets"
	#print "Does the euler characteristic hold?" + len(Pts) + len(Facets) == len(Edges)
	"""
	print ""
	print ""
