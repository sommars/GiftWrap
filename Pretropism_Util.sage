load("GiftWrap.sage")

def NormVector(V):
	Denom = 0
	for i in xrange(len(V)):
		Denom += V[i]^2
	Denom = Denom^(1/2)
	return [V[i]/Denom for i in xrange(len(V))]

#-------------------------------------------------------------------------------
def InnerProduct(V1,V2):
	InnerProd = 0
	for i in xrange(len(V1)):
		InnerProd += V1[i]*V2[i]
	return InnerProd

#-------------------------------------------------------------------------------
def FindInitialForm(Pts, V):
	InForm = []
	MinimalIP = 'Test'
	for Pt in Pts:
		IP = InnerProduct(V,Pt)
		if MinimalIP == 'Test':
			MinimalIP = IP
			InForm = [Pt]
		elif MinimalIP > IP:
			MinimalIP = IP
			InForm = [Pt]
		elif IP == MinimalIP:
			InForm.append(Pt)
	return InForm

#-------------------------------------------------------------------------------
def ComputeFacetPretropisms(AFaces, BFaces):
	#Currently is an inefficient implementation
	AFacetVectors = [(NormVector(AFaces[1][i].InnerNormals[0]),i) for i in xrange(len(AFaces[1]))]
	BFacetVectors = [NormVector(BFaces[1][i].InnerNormals[0]) for i in xrange(len(BFaces[1]))]

	FacetPretropisms = []
	NotFacetPretropisms = []
	for Vector in AFacetVectors:
		if Vector[0] in BFacetVectors:
			FacetPretropisms.append(Vector)
		else:
			NotFacetPretropisms.append(Vector)
	return FacetPretropisms, NotFacetPretropisms
"""
def ComputeFacetPretropisms2(AFaces, BFaces):
	AFacetVectors = [(NormVector(AFaces[1][i].InnerNormal[0]),i) for i in xrange(len(AFaces[1]))]
	BFacetVectors = [NormVector(BFaces[1][i].InnerNormal[0]) for i in xrange(len(BFaces[1]))]
	AFacetVectors.sort()
	print AFacetVectors
	FacetPretropisms = []
	NotFacetPretropisms = []
	i = 0
	for i in xrange(len(AFacetVectors)):
		Vector = AFacetVectors[i]
		while True:
			if len(BFacetVectors) == 0:
				return FacetPretropisms, NotFacetPretropisms + [AFacetVectors[j] for j in xrange(i,len(AFacetVectors))]
			if Vector[0] > BFacetVectors[0]:
				BFacetVectors.pop(0)
			else:
				break
		if Vector[0] < BFacetVectors[i]:
			NotFacetPretropisms.append(Vector)
		else:
			FacetPretropisms.append(Vector)
	return FacetPretropisms, NotFacetPretropisms
"""

#-------------------------------------------------------------------------------
def ConesDoIntersect(VList1, VList2):
	C1 = Cone(VList1)
	C2 = Cone(VList2)
	if len(C1.intersection(C2).rays()) == 0:
		return False
	else:
		return True

#-------------------------------------------------------------------------------
