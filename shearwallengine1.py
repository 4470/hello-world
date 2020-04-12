from __future__ import division
print " Please wait engine loads"
import ifcopenshell
import pandas
import numpy
import csv
import math

round_off = lambda x:int(x+0.5)
inp_file = raw_input("Enter the file without extension:")
ff = ifcopenshell.open(inp_file+".ifc")
xel = pandas.read_excel(r'c:\Users\narindersingh\Desktop\Narinder\masonary_sheet.xlsx')
walls = ff.by_type("ifcwall")
rel_material = ff.by_type("ifcrelassociatesmaterial")

no_of_walls = len(walls)

rein_in_wall_list = [[] for _ in range(no_of_walls)]

reins = ff.by_type("ifcreinforcingbar")  # All bundles of reinforcement

out_file = open(r"c:\Users\narindersingh\Desktop\Narinder\masonary.txt","w+")
disk = csv.writer(out_file)

disk.writerow(['................................................................... '])
disk.writerow(['        Report of Code Checks on Shear Wall as per IS:13920      '])
disk.writerow(['                Developed by GCE, GNDEC, Ludhiana                '])
disk.writerow(['...................................................................'])

print "File is compiled. Check report at result file generated"

for rein in reins:
    rods = rein.Representation.Representations[0].Items  # Rods in a list
    rod_mid = rods[round_off(len(rods)/2)]
    rod_mid_coord_i = rod_mid.Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates
    rod_mid_coord_f = rod_mid.Outer.CfsFaces[0].Bounds[0].Bound.Polygon[1].Coordinates

    rod_mid_coord_i_x = rod_mid_coord_i[0]
    rod_mid_coord_i_y = rod_mid_coord_i[1]
    rod_mid_coord_i_z = rod_mid_coord_i[2]
    
    rod_mid_coord_f_x = rod_mid_coord_f[0]
    rod_mid_coord_f_y = rod_mid_coord_f[1]
    rod_mid_coord_f_z = rod_mid_coord_f[2]

    for i in range(len(walls)):
        
        
        dir = walls[i].Representation.Representations[0].Items[0].Position.RefDirection
        x_ref_wall=round_off(dir[0][0])
        y_ref_wall=round_off(dir[0][1])
        z_ref_wall=round_off(dir[0][2])

        cart_point = walls[i].Representation.Representations[0].Items[0].Position.Location[0]
        x_cart_wall = cart_point[0]
        y_cart_wall = cart_point[1]
        z_cart_wall = cart_point[2]

        wall_thick = walls[i].Representation.Representations[0].Items[0].SweptArea.YDim

        if x_ref_wall == 0:
            lb=(x_cart_wall - wall_thick/2)
            ub=(x_cart_wall + wall_thick/2)
            if (lb<= rod_mid_coord_i_x <= ub) and (lb<= rod_mid_coord_f_x <= ub):
                rein_in_wall_list[i].append(rein)
                #print(rein.Name + " is in "+ walls[i].Name)
                
        elif y_ref_wall == 0:
            lb=(y_cart_wall - wall_thick/2)
            ub=(y_cart_wall + wall_thick/2)
            if (lb<= rod_mid_coord_i_y <=ub) and (lb<= rod_mid_coord_f_y <=ub):
                rein_in_wall_list[i].append(rein)
                #print(rein.Name + " is in "+ walls[i].Name)


for q, wall, array in zip(range(len(walls)), walls, rein_in_wall_list):
    for j in range(len(walls)):                                                                              #extract structural data
        if xel['walls_name'][j] == int(wall.Name):           
            load_induced = xel['load_induced'][j]
            sf_induced = xel['sf_induced'][j]
            bm_induced = xel['bm_induced'][j]
            load_seismic = xel['load_seismic'][j]
            sf_seismic = xel['sf_seismic'][j]
            bm_seismic = xel['bm_seismic'][j]

    Pa = ((0.8*load_induced)+(1.2*load_seismic))
    Pb = ((1.2*load_induced)+(1.2*load_seismic))
    Mu = (1.2*(bm_induced+bm_seismic))
    Vu = (1.2*(sf_induced+sf_seismic))

    checks_passed = []

    wall_len = wall.Representation.Representations[0].Items[0].SweptArea.XDim                               #extract wall dimentions
    wall_thk = wall.Representation.Representations[0].Items[0].SweptArea.YDim
    wall_ht = wall.Representation.Representations[0].Items[0].Depth

    for no in range(len(rel_material)):
        for components in rel_material[no]:
            for intels in components:
                if str(intels) == str(wall):                                                                #extract wall material
                    conc_mat = rel_material[no].RelatingMaterial.Name
    fck = int(conc_mat.split('M')[1])

    for no in range(len(rel_material)):
        for components in rel_material[no]:
            for intels in components:
                if str(intels) == str(array[0]):                                                           #extract reinforcement material
                    steel_mat = rel_material[no].RelatingMaterial.Name
    fy = int(steel_mat.split('Fe')[1])
    

    MOI = (((wall_thk*1000)*((wall_len*1000)**3))/12)
    G_area = ((wall_thk*1000)*(wall_len*1000))                      # dimentions in mm                 
    Z_mod = ((MOI)/((wall_len*1000)/(2)))

    fc_max = (((Pb*1000)/(G_area))+((Mu*1000000)/(Z_mod)))
    fc_min = (((Pb*1000)/(G_area))-((Mu*1000000)/(Z_mod)))
    fc_allow = (0.2*fck)

    disk.writerow([' Report for wall: '+ str(q + 1) ])

    
    # now check whether wall thickness is adequate
    if wall_thk >= 0.15:
        disk.writerow([' 1. Thickness of every part of wall is more than 150 mm: Check passed by '+ str((wall_thk*1000)-(150)) + 'mm'])
        checks_passed.append(str(1))
    else:
        disk.writerow([' 1. Thickness of every part of wall is more than 150 mm: Check failed by '+ str((150)-(wall_thk*1000)) + 'mm'])

    
    # also check no wall shoul be less than 150 mm
    if fc_allow >= fc_max:
        disk.writerow([' 2. Wall thickness is sufficient as per extreme fibre stresses '])
        checks_passed.append(str(2))
    else:
        disk.writerow([' 2. Wall thickness is not sufficient as per extreme fibre stresses '])

    
    eff_dep = (0.9*wall_len)
    
    tau_v = ((Vu*1000)/(wall_thk*1000*eff_dep*1000))                                                                     
    tau_v_lim = (0.25*(math.sqrt(fck)))

    # if wall thickness is more than 200 mm than provide reinforcements in 2 curtains
    if wall_thk >= .2 or tau_v >= tau_v_lim:
        if len(array) == 4:
            disk.writerow([' 3. Reinforcements should be provided in two curtains: Check passed '])
            checks_passed.append(str(3))
        else:
            disk.writerow([' 3. Reinforcements should be provided in two curtains: Check failed '])
    else:
        disk.writerow([' 3. Reinforcements provided in any curtains: Check passed '])
        checks_passed.append(str(3))
    

    for bundle in array:
        az = bundle.Representation.Representations[0].Items[0].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[2]
        bz = bundle.Representation.Representations[0].Items[0].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[1].Coordinates[2]
        if (abs((bz-az)*1000)) >= (bundle.BarLength-(0.05*(bundle.BarLength))):
            longi_bar = bundle
        else:
            trans_bar = bundle


    steel_bars = longi_bar.Representation.Representations[0].Items
    ax_longi = steel_bars[(len(steel_bars))//2].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[0]
    ay_longi = steel_bars[(len(steel_bars))//2].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[1]
    az_longi = steel_bars[(len(steel_bars))//2].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[2]
    bx_longi = steel_bars[((len(steel_bars))//2)+1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[0]
    by_longi = steel_bars[((len(steel_bars))//2)+1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[1]
    bz_longi = steel_bars[((len(steel_bars))//2)+1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[2]                 #spacing of centre longitudal bars
    centre_spacing = math.sqrt(((bx_longi - ax_longi)**2)+((by_longi - ay_longi)**2)+((bz_longi - az_longi)**2))

    dia_longi = longi_bar.NominalDiameter
    Pt_longi_oneline = (((0.785*((dia_longi)**2))/(centre_spacing*1000000*wall_thk))*100)

    # check if 2*Pt_longi_oneline > .25 ie. minimum steel check
    if (2*Pt_longi_oneline) >= .25:
        disk.writerow([' 4. Minimum longitudinal reinforcement: Check passed by '+ str(round(((((2*Pt_longi_oneline)-(0.25))/(.25))*100),2)) +'%' ])
        checks_passed.append(str(4))
    else:
        disk.writerow([' 4. Minimum longitudinal reinforcement: Check failed by '+ str(round(((((0.25)-(2*Pt_longi_oneline))/(.25))*100),2)) +'%' ])
        

    steel_bars_2 = trans_bar.Representation.Representations[0].Items
    ax_trans = steel_bars_2[(len(steel_bars))//2].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[0]
    ay_trans = steel_bars_2[(len(steel_bars))//2].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[1]
    az_trans = steel_bars_2[(len(steel_bars))//2].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[2]
    bx_trans = steel_bars_2[((len(steel_bars))//2)+1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[0]
    by_trans = steel_bars_2[((len(steel_bars))//2)+1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[1]
    bz_trans = steel_bars_2[((len(steel_bars))//2)+1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[2]
    trans_spacing = math.sqrt(((bx_trans - ax_trans)**2)+((by_trans - ay_trans)**2)+((bz_trans - az_trans)**2))

    dia_trans = trans_bar.NominalDiameter
    Pt_trans_oneline = (((0.785*((dia_trans)**2))/(trans_spacing*1000000*wall_thk))*100)

    # check if 2*Pt_transi_oneline > .25 ie. minimum steel check
    if (2*Pt_trans_oneline) >= .25:
        disk.writerow([' 5. Minimum transverse reinforcement: Check passed by '+ str(round(((((2*Pt_trans_oneline)-(0.25))/(.25))*100),2)) +'%' ])
        checks_passed.append(str(5))
    else:
        disk.writerow([' 5. Minimum transverse reinforcement: Check failed by '+ str(round(((((0.25)-(2*Pt_trans_oneline))/(.25))*100),2)) +'%' ])

    # Max. dia check
    if dia_longi <= (wall_thk*100):
        disk.writerow([' 6. Maximum diameter of reinforcement: Check passed by '+ str((wall_thk*100)-(dia_longi)) + 'mm' ])
        checks_passed.append(str(6))
    else:
        disk.writerow([' 6. Maximum diameter of reinforcement: Check failed by '+ str((dia_longi)-(wall_thk*100)) + 'mm' ])

    # Max spacing check
    if centre_spacing <= (wall_len/5) and centre_spacing <= (wall_thk*3) and centre_spacing <= (0.45):
        disk.writerow([' 7. Maximum spacing of verticle reinforcement: Check passed by '+ str(round(((450)-(centre_spacing*1000)),2)) + 'mm' ])
        checks_passed.append(str(7))
    else:
        disk.writerow([' 7. Maximum spacing of verticle reinforcement: Check failed by '+ str(round(((centre_spacing*1000)-(450)),2)) + 'mm' ])

    if trans_spacing <= (wall_len/5) and trans_spacing <= (wall_thk*3) and trans_spacing <= (0.45):
        disk.writerow([' 8. Maximum spacing of horizontal reinforcement: Check passed by '+ str(round(((450)-(trans_spacing*1000)),2)) + 'mm' ])
        checks_passed.append(str(8))
    else:
        disk.writerow([' 8. Maximum spacing of horizontal reinforcement: Check failed by '+ str(round(((trans_spacing*1000)-(450)),2)) + 'mm' ])

    
    Pt_longi = (Pt_longi_oneline*2)
    
    pt = [0.15, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00]

    Mi15 = [0.28, 0.35, 0.46, 0.54, 0.60, 0.64, 0.68, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71]
    Mi20 = [0.28, 0.36, 0.48, 0.56, 0.62, 0.66, 0.72, 0.75, 0.79, 0.81, 0.82, 0.82, 0.82]
    Mi25 = [0.29, 0.36, 0.49, 0.57, 0.64, 0.70, 0.74, 0.78, 0.82, 0.85, 0.88, 0.90, 0.92]
    Mi30 = [0.29, 0.37, 0.50, 0.59, 0.66, 0.71, 0.76, 0.80, 0.84, 0.88, 0.91, 0.94, 0.96]
    Mi35 = [0.29, 0.37, 0.50, 0.59, 0.67, 0.73, 0.78, 0.82, 0.86, 0.90, 0.93, 0.96, 0.99]
    Mi40 = [0.30, 0.38, 0.51, 0.60, 0.68, 0.74, 0.79, 0.84, 0.88, 0.92, 0.95, 0.98, 1.01]

    pt_cal = 1.27
    if Pt_longi <= 0.15:
        tau_c = 0.28
    else:
        if conc_mat == 'M15':
            tau_c = numpy.interp(Pt_longi,pt,Mi15)
        elif conc_mat == 'M20':
            tau_c = numpy.interp(Pt_longi,pt,Mi20)
        elif conc_mat == 'M25':
            tau_c = numpy.interp(Pt_longi,pt,Mi25)
        elif conc_mat == 'M30':
            tau_c = numpy.interp(Pt_longi,pt,Mi30)
        elif conc_mat == 'M35':
            tau_c = numpy.interp(Pt_longi,pt,Mi35)
        elif conc_mat == 'M40':
            tau_c = numpy.interp(Pt_longi,pt,Mi40)


    V_us = ((tau_v - tau_c)*(wall_thk*1000*1000))
    as_trans = (0.785*((dia_trans)**2))
    sv_req = ((0.87*fy*2*as_trans*1000)/(V_us))
    Vus_obt = ((0.87*fy*2*as_trans*1000)/(trans_spacing*1000))
    # chech shear ie. if sv_req < trans_spacing*1000 then safe
    if Vus_obt <= V_us:
        disk.writerow([' 9. Shear check failed by: '+ str(round((((V_us)-(Vus_obt))/1000),2)) + 'Kn' ])
    else:
        disk.writerow([' 9. Shear check passed by: '+ str(round((((Vus_obt)-(V_us))/1000),2)) + 'Kn' ])
        checks_passed.append(str(9))
    

    P_on_wall = (0.8*Pa)
    lamda = ((P_on_wall*1000)/(fck*wall_thk*1000*0.8*wall_len*1000))
    phi = ((0.87*fy*Pt_longi)/(100*fck))
    beta = ((0.87*fy)/700)
    lim_naxis = (0.0035/(0.0035+((0.87*fy)/200000)))
    naxis_try = ((phi+lamda)/((2*phi)+0.36))

    if naxis_try <= lim_naxis :
        naxis = naxis_try
        Muv_wall = ((fck*wall_thk*1000*((0.8*wall_len*1000)**2)*phi)*(((1+(lamda/phi))*(0.5-(0.416*naxis)))-((naxis**2)*(0.168+((beta**2)/3)))))

    elif naxis_try >= lim_naxis :
        alpa_1 = (0.36+(phi*(1-(0.5*beta)-(1/(2*beta)))))
        alpa_2 = (0.15+((phi/2)*(1-beta-((beta**2)/2)-(1/(3*beta)))))
        alpa_4 = ((phi/beta)-lamda)
        alpa_5 = (phi/(2*beta))
        naxis = ((-alpa_4+(math.sqrt((alpa_4**2)-(4*alpa_1*alpa_5))))/(2*alpa_1))
        Muv_wall = ((fck*wall_thk*1000*((0.8*wall_len*1000)**2))*((alpa_1*naxis)-(alpa_2*(naxis**2))-(alpa_3)-(lamda/2)))
        
    axuu = longi_bar.Representation.Representations[0].Items[0].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[0]
    ayuu = longi_bar.Representation.Representations[0].Items[0].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[1]
    azuu = longi_bar.Representation.Representations[0].Items[0].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[2]
    bxuu = longi_bar.Representation.Representations[0].Items[1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[0]
    byuu = longi_bar.Representation.Representations[0].Items[1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[1]
    bzuu = longi_bar.Representation.Representations[0].Items[1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[2]
    spacing_beg = math.sqrt(((bxuu-axuu)**2)+((byuu-ayuu)**2)+((bzuu-azuu)**2))

    lb_spacing = spacing_beg - (0.05*spacing_beg)
    ub_spacing = spacing_beg + (0.05*spacing_beg)

    momentBars = []
    todoreinbars = longi_bar.Representation.Representations[0].Items
    enespaces = len(todoreinbars) - 1
    for g in range(enespaces):
        ao_x = todoreinbars[g].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[0]
        ao_y = todoreinbars[g].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[1]
        ao_z = todoreinbars[g].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[2]
        bo_x = todoreinbars[g+1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[0]
        bo_y = todoreinbars[g+1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[1]
        bo_z = todoreinbars[g+1].Outer.CfsFaces[0].Bounds[0].Bound.Polygon[0].Coordinates[2]
        spacing_nik = math.sqrt(((bo_x-ao_x)**2)+((bo_y-ao_y)**2)+((bo_z-ao_z)**2))
        if lb_spacing <= spacing_nik <= ub_spacing:
            momentBars.append(todoreinbars[g])

    no_end_bars = (2+len(momentBars))
    as_endbar = (0.785*(dia_longi**2))
    Pu_boundry_ele = ((0.4*fck*((0.1*wall_len*1000*wall_thk*1000)-(no_end_bars*as_endbar)))+(0.67*fy*(no_end_bars*as_endbar)))
    bdy_mo_cap = (Pu_boundry_ele*0.9*wall_len*1000)
    Mu_cap_all = ( bdy_mo_cap + Muv_wall)
    
    # check moment capacity
    if Mu_cap_all >= (Mu*1000000):
        disk.writerow([' 10. Moment check passed by: '+ str(round(((Mu_cap_all-(Mu*1000000))/1000000),2))+' Kn-m' ])
        checks_passed.append(str(10))
    else:
        disk.writerow([' 10. Moment check failed by: '+ str(round((((Mu*1000000)-Mu_cap_all)/1000000),2))+' Kn-m' ])
    
    
    disk.writerow(['\n'])
    disk.writerow(['               Compliance result for shear wall: '+ str((len(checks_passed))*10) +'% checks are passed'])

    

    
    disk.writerow(['\n\n\n'])         
    


            


    
    

    
    
    
    


    
    
    
    




    print ("\n\n\n")


out_file.close()    
    
    
    
    


    


    
    


    


    

    


    


    
    

            
            
    
    
    
    
    

