#This will only work for n x (n+1)
#assumes the matrix is coming in as a list of lists, but already in Hermite Normal Form

#-------------------------------------------------------------------------------
def FirstNonZeroEntryOfRowAndItsIndex(Row):
	"""
	Finds the first value in a matrix that is not zero. The first return is the 
	boolean value RowIsZero, the second is the entry itself, and the third is the
	index where it is not zero
	"""
	for i in range(len(Row)):
		if Row[i] != 0:
			return False, Row[i], i
	return True, -1, -1

#-------------------------------------------------------------------------------
def RowIsZero(Row):
	"""
	Function that returns a boolean based on whether or not a row is all zeroes
	"""
	for Int in Row:
		if Int != 0:
			return False
	return True

#-------------------------------------------------------------------------------
def lcmm(x, y):
	if x > y:
		greater = x
	else:
		greater = y
	while(True):
		if((greater % x == 0) and (greater % y == 0)):
			lcm = greater
			break
			greater += 1
	return lcm

#-------------------------------------------------------------------------------
def gcd(a, b):
	"""Return greatest common divisor using Euclid's Algorithm."""
	while b:
		a, b = b, a % b
	return a

#-------------------------------------------------------------------------------
def lcm(a, b):
	"""Return lowest common multiple."""
	return a * b // gcd(a, b)

#-------------------------------------------------------------------------------
def DotProductPartialRowsStartingAtIndex(Row, Index):
	#print NormalVector, Row, Index
	Sum = 0
	for i in range(Index, len(Row)):
		Sum += NormalVector[i]*Row[i]
	return Sum

#-------------------------------------------------------------------------------
def InitializeNormalVectorForGivenIndex(Index):
	#print "asf", NormalVector, Index
	for i in range(Index, len(NormalVector)):
		if NormalVector[i] == 'n':
			NormalVector[i] = 1
	return

#-------------------------------------------------------------------------------
def MultiplyNormalVectorByValue(Value):
	for i in range(len(NormalVector)):
		if NormalVector[i] != 'n':
			NormalVector[i] = NormalVector[i]*Value
	return

#-------------------------------------------------------------------------------
def RefineNormal(Row, Index):
	if RowIsZero(Row):
		return
	InitializeNormalVectorForGivenIndex(Index+1)
	Sum = DotProductPartialRowsStartingAtIndex(Row, Index+1)
	if Sum%Row[Index] != 0:
		Templcm = lcm(Sum, Row[Index])
		MultiplyNormalVectorByValue(Templcm/Sum)
		Sum = Sum*Templcm
		NormalVector[Index] = -Templcm/Row[Index]
	else:
		NormalVector[Index] = -Sum/Row[Index]
	return

#-------------------------------------------------------------------------------
def CheckNormal(Rows):
	if RowIsZero(NormalVector):
		print "The normal vector is zero. This means something is broken"
		return
	NewVector = []
	for Row in Rows:
		TempInt = 0
		for i in range(len(Row)):
			TempInt += Row[i]*NormalVector[i]
		NewVector.append(TempInt)
	#print NewVector
	if not RowIsZero(NewVector):
		#print "It works"
		#print NormalVector
	#else:
		print "SOMETHING IS BROKEN"
		print "here is the matrix"
		Print(Rows)
		print "here is your normal"
		Print(NormalVector)
	return

#-------------------------------------------------------------------------------
def Print(Rows):
	for Row in Rows:
		print Row
	print ""
	return

#-------------------------------------------------------------------------------
def ColumnOfZeroGivesNormalVectorOne(Rows):
	for i in range(len(Rows[0])):
		for j in range(len(Rows)):
			if Rows[j][i] != 0:
				break
			elif j == len(Rows) - 1:
				NormalVector[i] = 1
	return

#-------------------------------------------------------------------------------
def ReduceRow(Row):
	RowIsZero, RowEntry, Index = FirstNonZeroEntryOfRowAndItsIndex(Row)
	if RowIsZero == True:
		return Row
	if RowEntry == 1:
		return Row
	CommonFactor = RowEntry
	NumberOfNonZeroTerms = 1
	for i in range(Index + 1, len(Row)):
		CommonFactor = egcd(CommonFactor, Row[i])[0]
		if CommonFactor == 1:
			return Row
		if Row[i] != 0:
			NumberOfNonZeroTerms += 1
	if NumberOfNonZeroTerms < 2:
		return Row
	ReducedRow = []
	for Int in Row:
		ReducedRow.append(Int/CommonFactor)
	return ReducedRow

#-------------------------------------------------------------------------------
def FindNormal(Rows):
	global NormalVector
	NormalVector = []
	for k in range(len(Rows[0])):
		NormalVector.append('n')

	for i in range(len(Rows)):
		#note to self change sentinel value of index to -1
		Index = len(Rows) - i - 1
		RowIsZero, RowEntry, RowIndex = FirstNonZeroEntryOfRowAndItsIndex(Rows[Index])
		if RowIndex != - 1:
			RefineNormal(Rows[Index], RowIndex)
	ColumnOfZeroGivesNormalVectorOne(Rows)
	CheckNormal(Rows)
	return NormalVector
