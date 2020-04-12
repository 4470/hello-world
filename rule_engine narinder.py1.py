import ifcopenshell
ifc_file = ifcopenshell.open("C:\Users\narindersingh\Desktop\Narinder\sdfghj.ifc")
walls = ifc_file.by_type("ifcwalls")
round_off = lambda x:int(x+0.5)
print('walls')
i = 1
cond1 = int(input("1: "))
cond2 = int(input("2: "))
cond3 = int(input("3: "))
G1 = 0
effective_height = 0
ff = ifcopenshell.open("C:\Users\narindersingh\Desktop\Narinder\a.ifc")
xel = pandas.read_excel(r'c:\Users\narindersingh\Desktop\Narinder\masonary_sheet.xlsx')
def calc_eff_height(cond3):
	if(cond3 == 1): G1 = 0.75*B4
	elif(cond3 == 2): G1 = 0.85*B4
	elif(cond3 == 3): G1 = 1*B4   
	elif(cond3 == 4): G1 = 1.5*B4
	elif(cond3 == 5): G1 = 0.75*B4 + 0.25*B6          
	elif(cond3 == 6 and B6 < 0.5*B4): G1 = B4
	else: G1 = 2*B4 
	if(cond3 == 5): effective_height = G1*B4
	elif(cond3 == 6): effective_height = G1*B4
	else: effective_height = G1*B4
	return effective_height
def calc_eff_length(cond1):
	if(cond1 == 1): Effective_length = 0.8*A2
	elif(cond1 == 2): Effective_length = 0.9*A2
	elif(cond1 == 3): Effective_length = A2
	elif(cond1 == 4): Effective_length = 1.5*A2
	elif(cond1 == 5): Effective_length = 2*A2
	else: Effective_length = A2

	return Effective_length


for wall in walls:
	geometry = wall.Representation.Representations[0].Items[0] # An IfcExtrudedAreaSolid in this example
	l = geometry.SweptArea.XDim
	b = geometry.SweptArea.YDim
	h = geometry.Depth

#	if wall.hasOpenings():
#	print("Has Openings")
for wall in walls:
  for definition in wall.IsDefinedBy:
      property_set = definition.RelatingPropertyDefinition
      print(property_set) # Might return Pset_WallCommon

for i in range(len(walls)):	
	cart_point = walls[i].Representation.Representations[0].Items[0].Position.Location[0]
	print(walls[i].GlobalId , cart_point) 
	heff = calc_eff_height(cond3)
	print("Heff: " + str(heff))

	leff = calc_eff_length(cond1)
	print(l)
	print("Leff: " + str(leff))

	teff = b	
	
	slenderness_ratio1 = heff/teff
	slenderness_ratio2 = leff/teff
	print(l,b,h)
	for definition in wall.IsDefinedBy:
		print('hello')
		property_set = definition.RelatingPropertyDefinition
		print(property_set + "123456") # Might return Pset_WallCommon

	print(wall.IsDefinedBy)

    
    for i in range(len(walls)):
	dir = walls[i].Representation.Representations[0].Items[0].Position.RefDirection
	xrw=round_off(dir[0][0])
	yrw=round_off(dir[0][1])
	zrw=round_off(dir[0][2])
	print(xrw)
	print(yrw)
	print(zrw)
	print()


	for no in range(len(rel_material)):
	        for components in rel_material[no]:
        	    for intels in components:
                	if str(intels) == str(wall):                                          #extract wall material
                    		conc_mat = rel_material[no].RelatingMaterial.Name

	print("Wall" + str(i) + ": " + wall.GlobalId)
	print("Length: "+ str(l))
	print("Thickness: "+ str(b))
	print("Height: "+ str(h))
	print(conc_mat)
	print()
	i = i+1

i = 0

num = int(input("Enter wall no (from 0 to 41): "))
if(walls[num].HasAssociations):
	print(walls[num].HasAssociations)
	print()

else:
	print("Property dne")


if(walls[num].ContainedInStructure):
	print(walls[num].ContainedInStructure)
	print()

else:
	print("Property dne")

if(walls[num].HasOpenings):
	print(walls[num].HasOpenings)
	print()

else:
	print("Property dne")

if(walls[num].ObjectPlacement):
	print(walls[num].ObjectPlacement)
	print()

else:
	print("Property dne")

if(walls[num].Representation):
	print(walls[num].Representation)
	print()



while(i<len(walls)):
	print(walls[i].HasAssociations)
	print(walls[i].ContainedInStructure)
	print(walls[i].HasOpenings)
	print(walls[i].ObjectPlacement)
	print(walls[i].Representation)
	print()


while(i<len(walls)): 
	j = 0
	while(j<len(walls[i].ContainedInStructure[0][4])):
		print(walls[i].ContainedInStructure[0][4][j])
		j= j + 1
	print()
	i = i + 1





	
