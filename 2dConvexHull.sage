def FindNextPt(Pts, Pt1):
	def DoTurn(Pt1, Pt2, Pt3):
	 	return cmp((Pt2[X] - Pt1[X])*(Pt3[Y] - Pt1[Y]) - (Pt3[X] - Pt1[X])*(Pt2[Y] - Pt1[Y]), 0)
	def Distance(p, q):
		return (q[X] - p[X])^2 + (q[Y] - p[Y])^2
	Pt2 = Pt1
	for Pt3 in Pts:
		Turn = DoTurn(Pt1, Pt2, Pt3)
		if Turn == -1 or Turn == 0 and Distance(Pt1, Pt3) > Distance(Pt1, Pt2):
			Pt2 = Pt3
	return Pt2

#-------------------------------------------------------------------------------
def FindDistinctCoords(Pts):
	Coords = []
	for i in xrange(1,len(Pts)):
		for j in xrange(len(Pts[i])):
			Coord = Pts[i][j]
			if (Coord != Pts[0][j]) and (j not in Coords):
					Coords.append(j)
			if len(Coords) == 2:
				return Coords[0], Coords[1]
	return

#-------------------------------------------------------------------------------
def FindExtraPoints(Pts, Hull):
	ExtraPts = []
	for Pt in Pts:
		if Pt not in Hull:
			ExtraPts.append(Pt)
	return ExtraPts

#-------------------------------------------------------------------------------
def ConvexHull2d(Pts):
	Hull = [min(Pts)]
	global X
	global Y
	X, Y = 	FindDistinctCoords(Pts)
	for Pt in Hull:
		Pt2 = FindNextPt(Pts, Pt)
		if Pt2 != Hull[0]:
			Hull.append(Pt2)
	return Hull, FindExtraPoints(Pts,Hull)
