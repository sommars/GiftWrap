load("GiftUtil.sage")
# The code in this unit is the code that will have to change for going between 3 and n dimensions.
def GiftWrap(Pts):
	Facets = []
	if not PtsAreValid(Pts):
		return "The input set of points is not valid."
	Pts = RemoveDups(Pts)
	global Barycenter
	Barycenter = FindBarycenter(Pts)
	if CheckNormalFormDim(GetHNF(Pts)) != 2:
		print "This set of point is not full dimensional"
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
		FacetPts = Edge + FindNewFacetPtsFromEdge(Pts, Edge, Facet.InnerNormal, Facet.Vertices)

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

	# Euler Characteristic check
	EdgeCount = 0
	for Facet in Facets:
		EdgeCount += len(Facet.Vertices)
	# We don't want to double count edges
	EdgeCount = EdgeCount/2
	if len(Pts) - EdgeCount + len(Facets) != 2:
		print "Euler Characteristic failed!"
		raw_input()
	return Facets

#-------------------------------------------------------------------------------
def FindInitialFacet(Pts):
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

	if len(FirstPts) == 1:
		#Start with one point
		FirstPts = FirstPts + FindNewFacetPtsFromSinglePt(Pts, FirstPts)

		# Now we have either an edge or a facet. Need to check which via a dimension calculation
		if CheckNormalFormDim(GetHNF(FirstPts)) == 1:
			return MakeFacet(FirstPts, Pts)

		# Here I manufacture a point to create a fake facet. This allows me to
		# have a valid pair of normals, though they aren't ultimately meaningful
		# to the final convex hull
		KnownFacetPts = list(FirstPts)
		NewPoint = [FirstPts[0][0],FirstPts[0][1]+1,FirstPts[0][2]]
		if not NewPoint in KnownFacetPts:
			KnownFacetPts.append(NewPoint)
		InnerNormal = GetNormalFromHNF(GetHNF(KnownFacetPts))
		if not NormalPointsTowardsPt(InnerNormal, Barycenter, FirstPts[0]):
			InnerNormal = [-InnerNormal[i] for i in xrange(len(InnerNormal))]
		# Not sure why the integer calls below are necessary. Need to look into that
		for i in xrange(len(InnerNormal)):
			InnerNormal[i] = Integer(InnerNormal[i])

		FirstPts = FirstPts + FindNewFacetPtsFromEdge(Pts, FirstPts, InnerNormal, KnownFacetPts)
		return MakeFacet(FirstPts, Pts)		
	else:
		HNF = GetHNF(FirstPts)
		Dim = CheckNormalFormDim(HNF)
		# This is the case in which I immediately find a facet in vertical position
		if Dim == 1:
			return MakeFacet(FirstPts, Pts)

		# I have found an edge. I can do roughly the same thing I do in the general
		# situation. I have the same normal to the hyperplane [-1,0,0]. 
		# I must use a specific normal to the point set
		elif Dim == 0:
			Normal = [-1,0,0]
			
			# This is a way to manufacture a third point on the facet.
			KnownFacetPts = list(FirstPts)
			ExtraPt = [0 for i in xrange(len(FirstPts[0]))]
			ExtraPt[0] = FirstPts[0][0]
			for Pt in KnownFacetPts:
				for i in xrange(1,len(Pt)):
					ExtraPt[i] += abs(Pt[i])
			KnownFacetPts.append(ExtraPt)
			
			FirstPts = FirstPts + FindNewFacetPtsFromEdge(Pts, FirstPts, [-1,0,0], KnownFacetPts)
			return MakeFacet(FirstPts, Pts)
		else:
			print "Internal error, unexpected dimension = ", Dim
			raw_input()
			return
	return MakeFacet(FirstPts, Pts)

#-------------------------------------------------------------------------------
def MakeFacet(FacetPts, Pts):
	# We have all of the points in their general position. We need to bring them
	# into vertical position so that the convexhull algorithm can work.
	UCT = GetUCTAndNormal(GetNormalFromHNF(GetHNF(FacetPts)))[0]
	ShiftedFacetPts = TransformPts(FacetPts, UCT)
	Hull, PointsToRemove = ConvexHull2d(ShiftedFacetPts)
	
	# We want to make sure the hyperplane supporting the facet is to one side
	# of all of the points.
	CheckAllPtsLieOnOthersideOfFacetHyperplane(copy(Pts), copy(ShiftedFacetPts), UCT)

	Hull = TransformPts(Hull, matrix(UCT^-1, ZZ))
	
	#Hate doing this. Not sure why it's necessary, but this is a temporary solution
	for i in xrange(len(Hull)):
		for j in xrange(len(Hull[i])):
			Hull[i][j] = Integer(Hull[i][j])

	PointsToRemove = TransformPts(PointsToRemove, matrix(UCT^-1, ZZ))
	Facet = Face()

	Normal = GetNormalFromHNF(GetHNF(Hull))
	# Need to make the normal an inner normal
	if not NormalPointsTowardsPt(Normal, Barycenter, Hull[0]):
		Normal = [-Normal[i] for i in xrange(len(Normal))]
	Facet.InnerNormal = Normal
	Edges = [[Hull[len(Hull) - 1], Hull[0]]]
	for i in xrange(len(Hull) - 1):
		Edges.append([Hull[i],Hull[i+1]])
	Facet.Vertices = Hull
	Facet.Children = Hull
	Dim = CheckNormalFormDim(GetHNF(Hull))
	if Dim != 1:
		print "Internal error, 2d convex hull is in the wrong dimension."
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

def MakeRandomPointSet():
	Pts = []
	for i in xrange(randint(4,100)):
		Pts.append([])
		for j in xrange(3):
			Pts[i].append(Integer(randint(0,100)))
	return Pts

Tests = []
for i in xrange(5000):
	Tests.append((MakeRandomPointSet(), 'Random'))

#Tests = [(TestCube,'Cube'), (Tetrahedron,'Tetrahedron'),(TippedOverHouse,'TippedOverHouse'),(TetrahedronExtraPts,'TetrahedronExtraPts'),(BigCube,'BigCube')]

for Test in Tests:
	Pts = Test[0]
	print Test[1]
	print Pts
	GiftWrap(Pts)
	print "Done with " + Test[1]
	print ""
	print ""
