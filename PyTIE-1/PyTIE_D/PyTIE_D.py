class topological_indices_Degree:
    def topoexpressd():
        try:
            print("Provide your data to calculate Degree-based Topological Descriptors:")
            import numpy
            import math
            import sympy as sp
            x=sp.Symbol("x")
            y=sp.Symbol("y")
            from numpy.polynomial import Polynomial as P
            a1=[]
            a11=[]
            a111=[]
            b1=[]
            b3=[]
            b11=[]
            b111=[]
            c1=[]
            c11=[]
            c111=[]
            d111=[]
            f1=[]
            f2=[]
            f3=[]
            f4=[]
            f5=[]
            f6=[]
            f7=[]
            f8=[]
            f9=[]
            f10=[]
            f11=[]
            f12=[]
            f13=[]
            f14=[]
            f15=[]
            f16=[]
            f17=[]
            f18=[]
            f19=[]
            f20=[]
            f21=[]
            f22=[]
            f23=[]
            f24=[]
            f25=[]
            f1!=[0]
            f2!=[0]
            f3!=[0]
            f4!=[0]
            f5!=[0]
            f6!=[0]
            f7!=[0]
            f8!=[0]
            f9!=[0]
            f10!=[0]
            f11!=[0]
            f12!=[0]
            f13!=[0]
            f14!=[0]
            f15!=[0]
            f16!=[0]
            f17!=[0]
            f18!=[0]
            f19!=[0]
            f20!=[0]
            f21!=[0]
            f22!=[0]
            f23!=[0]
            f24!=[0]
            f25!=[0]
            h1=[]
            h2=[]
            h3=[]
            r1=[]
            #ranges of edge partition
            m=int(input("Specify the minimum degree of a given graph (G): ")) 
            n=int(input("Specify the maximum degree of a given graph (G): "))
            for i in range(m, n+1):
                for j in range(m, n+1):
                    for r in range(1, 7):
                        if(i>=j):
                            edge_part=int(input("Specify the values for partitioning edges in the graph of "+str(r)+" dimension of "+str((j, i))+": "))
                            a1.append(edge_part)
                            a11.append(edge_part)
                            a111.append(edge_part)
            #coding of expression of partition
            h=int(input("Specify the quantity of edge partitions generated above: "))
            for p in range(1, h+1):
                for i in range(((6*p)-5),6*p):
                    u0=(a1[i]-a1[i-1])
                    b1.append(u0)
                for i in range(((6*p)-5), 6*p):
                    u1=(a11[i]-a11[i-1])
                    b11.append(u1)
                for i in range(((5*p)-4),5*p):
                    u2=(b11[i]-b11[i-1])
                    c11.append(u2)
                for i in range(((6*p)-5),6*p):
                    u3=(a111[i]-a111[i-1])
                    b111.append(u3)
                for i in range(((5*p)-4), 5*p):
                    u4=(b111[i]-b111[i-1])
                    c111.append(u4)
                for i in range(((4*p)-3), 4*p):
                    u5=(c111[i]-c111[i-1])
                    d111.append(u5)
            print("Edge set partitions of a graph (G):")
            for q in range (1, h+1):
                if(a1[(6*q)-6]==a1[(6*q)-5]==a1[(6*q)-4]==a1[(6*q)-3]==a1[(6*q)-2]==a1[(6*q)-1]>=0):
                    p0=a1[(6*q)-6]
                    print(p0)
                    for i in range(m, n+1):
                        for j in range(m, n+1):
                            if(i>=j):
                                #First Zagreb index
                                x1=p0*((j+i))
                                f1.append(x1)
                                #Second Zagreb index
                                x2=p0*((j*i))
                                f2.append(x2)
                                #Hyper Zagreb index
                                x3=p0*((j+i)**2)
                                f3.append(x3)
                                #Third Zagreb index
                                x4=p0*(abs(j-i)) 
                                f4.append(x4)
                                #Reduced Zagreb index
                                x5=p0*((j-1)*(i-1))
                                f5.append(x5)
                                #Second modified Zagreb index
                                x6=p0*(1/(j*i)) 
                                f6.append(x6)
                                #Randic index
                                x7=p0*(1/math.sqrt(j*i)) 
                                f7.append(x7)
                                #Reciprocal Randic index
                                x8=p0*(math.sqrt(j*i))
                                f8.append(x8)
                                #Reduced reciprocal Randic index
                                x9=p0*(math.sqrt((j-1)*(i-1)))
                                f9.append(x9)
                                #General Randic index
                                x10=p0*((j*i)**(-1)) 
                                f10.append(x10)
                                #Atom bond connectivity index
                                x11=p0*(math.sqrt((j+i-2)/(j*i))) 
                                f11.append(x11)
                                #Geometric Arithematic index
                                x12=p0*(2*(math.sqrt(j*i)/(j+i)))
                                f12.append(x12)
                                #Harmonic index
                                x13=p0*(2/(j+i))
                                f13.append(x13)
                                #Sum-connectivity index
                                x14=p0*(1/(math.sqrt(j+i)))
                                f14.append(x14)
                                #Inverse sum index
                                x15=p0*((j*i)/(j+i))
                                f15.append(x15)
                                #Alberston index
                                x16=p0*(abs(j-i))
                                f16.append(x16)
                                #Symmetric division index
                                x17=p0*(((j**2)+(i**2))/(j*i)) 
                                f17.append(x17)
                                #Forgotten index
                                x18=p0*((j**2)+(i**2))
                                f18.append(x18)
                                #Sombor index
                                x19=p0*(math.sqrt((j**2)+(i**2)))
                                f19.append(x19)
                                #bi-Zagreb index
                                x20=p0*(j+i+(j*i))
                                f20.append(x20)
                                #Tri-Zagreb index
                                x21=p0*((j**2)+(i**2)+(j*i))
                                f21.append(x21)
                                #Geometric Harmonic index
                                x22=p0*((math.sqrt(j*i)*(j+i))/2)
                                f22.append(x22)
                                #Geometric Bi-Zagreb index
                                x23=p0*(math.sqrt(j*i)/(j+i+j*i))
                                f23.append(x23)
                                #Geometric Tri-Zagreb index
                                x24=p0*(math.sqrt(j*i)/((j**2)+(i**2)+j*i))
                                f24.append(x24)
                                #Harmonic Geometric index
                                x25=p0*(2/((math.sqrt(j*i))*(j+i)))
                                f25.append(x25)
                elif(b1[(5*q)-5]==b1[(5*q)-4]!=0):
                    A1=b1[(5*q)-5]
                    A2=a1[(6*q)-6]
                    A3=A2-A1
                    p1=A3+A1*x
                    print(p1)
                    for i in range(m, n+1):
                        for j in range(m, n+1):
                            if(i>=j):
                                #First Zagreb index
                                x1=p1*((j+i))
                                f1.append(x1)
                                #Second Zagreb index
                                x2=p1*((j*i))
                                f2.append(x2)
                                #Hyper Zagreb index
                                x3=p1*((j+i)**2)
                                f3.append(x3)
                                #Third Zagreb index
                                x4=p1*(abs(j-i)) 
                                f4.append(x4)
                                #Reduced Zagreb index
                                x5=p1*((j-1)*(i-1))
                                f5.append(x5)
                                #Second modified Zagreb index
                                x6=p1*(1/(j*i)) 
                                f6.append(x6)
                                #Randic index
                                x7=p1*(1/math.sqrt(j*i)) 
                                f7.append(x7)
                                #Reciprocal Randic index
                                x8=p1*(math.sqrt(j*i))
                                f8.append(x8)
                                #Reduced reciprocal Randic index
                                x9=p1*(math.sqrt((j-1)*(i-1)))
                                f9.append(x9)
                                #General Randic index
                                x10=p1*((j*i)**(-1)) 
                                f10.append(x10)
                                #Atom bond connectivity index
                                x11=p1*(math.sqrt((j+i-2)/(j*i))) 
                                f11.append(x11)
                                #Geometric Arithematic index
                                x12=p1*(2*(math.sqrt(j*i)/(j+i)))
                                f12.append(x12)
                                #Harmonic index
                                x13=p1*(2/(j+i))
                                f13.append(x13)
                                #Sum-connectivity index
                                x14=p1*(1/(math.sqrt(j+i)))
                                f14.append(x14)
                                #Inverse sum index
                                x15=p1*((j*i)/(j+i))
                                f15.append(x15)
                                #Alberston index
                                x16=p1*(abs(j-i))
                                f16.append(x16)
                                #Symmetric division index
                                x17=p1*(((j**2)+(i**2))/(j*i)) 
                                f17.append(x17)
                                #Forgotten index
                                x18=p1*((j**2)+(i**2))
                                f18.append(x18)
                                #Sombor index
                                x19=p1*(math.sqrt((j**2)+(i**2)))
                                f19.append(x19)
                                #bi-Zagreb index
                                x20=p1*(j+i+(j*i))
                                f20.append(x20)
                                #Tri-Zagreb index
                                x21=p1*((j**2)+(i**2)+(j*i))
                                f21.append(x21)
                                #Geometric Harmonic index
                                x22=p1*((math.sqrt(j*i)*(j+i))/2)
                                f22.append(x22)
                                #Geometric Bi-Zagreb index
                                x23=p1*(math.sqrt(j*i)/(j+i+j*i))
                                f23.append(x23)
                                #Geometric Tri-Zagreb index
                                x24=p1*(math.sqrt(j*i)/((j**2)+(i**2)+j*i))
                                f24.append(x24)
                                #Harmonic Geometric
                                x25=p1*(2/((math.sqrt(j*i))*(j+i)))
                                f25.append(x25)
                elif(c11[(4*q)-4]==c11[(4*q)-3]!=0):
                    B1=(1/2)*(c11[(4*q)-4])
                    B2=b11[(5*q)-5]-3*B1
                    B3=a11[(6*q)-6]-(B1+B2)
                    p2=B1*x**2+B2*x+B3
                    print(p2)
                    for i in range(m, n+1):
                        for j in range(m, n+1):
                            if(i>=j):
                                #First Zagreb index
                                x1=p2*((j+i))
                                f1.append(x1)
                                #Second Zagreb index
                                x2=p2*((j*i))
                                f2.append(x2)
                                #Hyper Zagreb index
                                x3=p2*((j+i)**2)
                                f3.append(x3)
                                #Third Zagreb index
                                x4=p2*(abs(j-i)) 
                                f4.append(x4)
                                #Reduced Zagreb index
                                x5=p2*((j-1)*(i-1))
                                f5.append(x5)
                                #Second modified Zagreb index
                                x6=p2*(1/(j*i)) 
                                f6.append(x6)
                                #Randic index
                                x7=p2*(1/math.sqrt(j*i)) 
                                f7.append(x7)
                                #Reciprocal Randic index
                                x8=p2*(math.sqrt(j*i))
                                f8.append(x8)
                                #Reduced reciprocal Randic index
                                x9=p2*(math.sqrt((j-1)*(i-1)))
                                f9.append(x9)
                                #General Randic index
                                x10=p2*((j*i)**(-1)) 
                                f10.append(x10)
                                #Atom bond connectivity index
                                x11=p2*(math.sqrt((j+i-2)/(j*i))) 
                                f11.append(x11)
                                #Geometric Arithematic index
                                x12=p2*(2*(math.sqrt(j*i)/(j+i)))
                                f12.append(x12)
                                #Harmonic index
                                x13=p2*(2/(j+i))
                                f13.append(x13)
                                #Sum-connectivity index
                                x14=p2*(1/(math.sqrt(j+i)))
                                f14.append(x14)
                                #Inverse sum index
                                x15=p2*((j*i)/(j+i))
                                f15.append(x15)
                                #Alberston index
                                x16=p2*(abs(j-i))
                                f16.append(x16)
                                #Symmetric division index
                                x17=p2*(((j**2)+(i**2))/(j*i)) 
                                f17.append(x17)
                                #Forgotten index
                                x18=p2*((j**2)+(i**2))
                                f18.append(x18)
                                #Sombor index
                                x19=p2*(math.sqrt((j**2)+(i**2)))
                                f19.append(x19)
                                #bi-Zagreb index
                                x20=p2*(j+i+(j*i))
                                f20.append(x20)
                                #Tri-Zagreb index
                                x21=p2*((j**2)+(i**2)+(j*i))
                                f21.append(x21)
                                #Geometric Harmonic index
                                x22=p2*((math.sqrt(j*i)*(j+i))/2)
                                f22.append(x22)
                                #Geometric Bi-Zagreb index
                                x23=p2*(math.sqrt(j*i)/(j+i+j*i))
                                f23.append(x23)
                                #Geometric Tri-Zagreb index
                                x24=p2*(math.sqrt(j*i)/((j**2)+(i**2)+j*i))
                                f24.append(x24)
                                #Harmonic Geometric
                                x25=p2*(2/((math.sqrt(j*i))*(j+i)))
                                f25.append(x25)
                elif(d111[(3*q)-3]==d111[(3*q)-2]!=0):
                    C1=(1/6)*(d111[3*q-3])
                    C2=(1/2)*(c111[(4*q)-4]-12*C1) 
                    C3=b111[(5*q)-5]-(3*C2+7*C1)
                    C4=a111[(6*q)-6]-(C1+C2+C3)
                    p3=C1*x**3+C2*x**2+C3*x+C4
                    print(p3)
                    for i in range(m, n+1):
                        for j in range(m, n+1):
                            if(i>=j):
                                #First Zagreb index
                                x1=p3*((j+i))
                                f1.append(x1)
                                #Second Zagreb index
                                x2=p3*((j*i))
                                f2.append(x2)
                                #Hyper Zagreb index
                                x3=p3*((j+i)**2)
                                f3.append(x3)
                                #Third Zagreb index
                                x4=p3*(abs(j-i)) 
                                f4.append(x4)
                                #Reduced Zagreb index
                                x5=p3*((j-1)*(i-1))
                                f5.append(x5)
                                #Second modified Zagreb index
                                x6=p3*(1/(j*i)) 
                                f6.append(x6)
                                #Randic index
                                x7=p3*(1/math.sqrt(j*i)) 
                                f7.append(x7)
                                #Reciprocal Randic index
                                x8=p3*(math.sqrt(j*i))
                                f8.append(x8)
                                #Reduced reciprocal Randic index
                                x9=p3*(math.sqrt((j-1)*(i-1)))
                                f9.append(x9)
                                #General Randic index
                                x10=p3*((j*i)**(-1)) 
                                f10.append(x10)
                                #Atom bond connectivity index
                                x11=p3*(math.sqrt((j+i-2)/(j*i))) 
                                f11.append(x11)
                                #Geometric Arithematic index
                                x12=p3*(2*(math.sqrt(j*i)/(j+i)))
                                f12.append(x12)
                                #Harmonic index
                                x13=p3*(2/(j+i))
                                f13.append(x13)
                                #Sum-connectivity index
                                x14=p3*(1/(math.sqrt(j+i)))
                                f14.append(x14)
                                #Inverse sum index
                                x15=p3*((j*i)/(j+i))
                                f15.append(x15)
                                #Alberston index
                                x16=p3*(abs(j-i))
                                f16.append(x16)
                                #Symmetric division index
                                x17=p3*(((j**2)+(i**2))/(j*i)) 
                                f17.append(x17)
                                #Forgotten index
                                x18=p3*((j**2)+(i**2))
                                f18.append(x18)
                                #Sombor index
                                x19=p3*(math.sqrt((j**2)+(i**2)))
                                f19.append(x19)
                                #bi-Zagreb index
                                x20=p3*(j+i+(j*i))
                                f20.append(x20)
                                #Tri-Zagreb index
                                x21=p3*((j**2)+(i**2)+(j*i))
                                f21.append(x21)
                                #Geometric Harmonic index
                                x22=p3*((math.sqrt(j*i)*(j+i))/2)
                                f22.append(x22)
                                #Geometric Bi-Zagreb index
                                x23=p3*(math.sqrt(j*i)/(j+i+j*i))
                                f23.append(x23)
                                #Geometric Tri-Zagreb index
                                x24=p3*(math.sqrt(j*i)/((j**2)+(i**2)+j*i))
                                f24.append(x24)
                                #Harmonic Geometric
                                x25=p3*(2/((math.sqrt(j*i))*(j+i)))
                                f25.append(x25)
            n_n_v=int(input("Enter the number of numerical values needed?: ")) #n_n_v=number of numerical values
            print("Python computational expressions for the calculation of topological indices: ")
            #COMPUTATION PART
            #First Zagreb computaion
            def first_Zagreb():
                sum1=0
                l1=[]
                l_n_v1=[]
                for q in range(1, h+1):
                    D1=f1[((h+1)*q)-(h+1)]
                    l1.append(D1)
                for i in l1:
                    sum1=sum1+i
                print("First Zagreb index =",(sum1))
                for i in range(1,n_n_v+1):
                    n_v=sum1.subs(x,i) #n_v-numerical values
                    l_n_v1.append(n_v) #l_n_v1-1st list of numerical values
                print('N.V of First Zagreb index = ', l_n_v1)
            first_Zagreb()
            #Second Zagreb computation
            def second_Zagreb():       
                sum2=0
                l2=[]
                l_n_v2=[]
                for q in range(1, h+1):
                    D2=f2[((h+1)*q)-(h+1)]
                    l2.append(D2)
                for i in l2:
                    sum2=sum2+i
                print("Second Zagreb index =",(sum2))
                for i in range(1,n_n_v+1):
                    n_v=sum2.subs(x,i) #n_v-numerical values
                    l_n_v2.append(n_v) #l_n_v2-2nd list of numerical values
                print('N.V of Second Zagreb index = ', l_n_v2)
            second_Zagreb()
            #Hyper Zagreb computation
            def hy_Zagreb(): 
                sum3=0
                l3=[]
                l_n_v3=[]
                for q in range(1, h+1):
                    D3=f3[((h+1)*q)-(h+1)]
                    l3.append(D3)
                for i in l3:
                    sum3=sum3+i
                print("Hyper Zagreb index =",(sum3))
                for i in range(1,n_n_v+1):
                    n_v=sum3.subs(x,i) #n_v-numerical values
                    l_n_v3.append(n_v) #l_n_v3-3rd list of numerical values
                print('N.V of Hyper Zagreb index = ', l_n_v3)
            hy_Zagreb() 
            #Third Zagreb index
            def third_Zagreb():     
                sum4=0
                l4=[]
                l_n_v4=[] #l_n_v4-4th list of numerical values
                for q in range(1, h+1):
                    D4=f4[((h+1)*q)-(h+1)]
                    l4.append(D4)
                for i in l4:
                    sum4=sum4+i
                print("Third Zagreb index =",(sum4))
                for i in range(1,n_n_v+1):
                    n_v=sum4.subs(x,i) #n_v-numerical values
                    l_n_v4.append(n_v) 
                print('N.V of Third Zagreb index = ', l_n_v4)
            third_Zagreb()
            #Reduced Zagreb index
            def redu_Zagreb():     
                sum5=0
                l5=[]
                l_n_v5=[] #l_n_v4-4th list of numerical values
                for q in range(1, h+1):
                    D5=f5[((h+1)*q)-(h+1)]
                    l5.append(D5)
                for i in l5:
                    sum5=sum5+i
                print("Reduced Zagreb index =",(sum5))
                for i in range(1,n_n_v+1):
                    n_v=sum5.subs(x,i) #n_v-numerical values
                    l_n_v5.append(n_v) 
                print('N.V of Reduced Zagreb index = ', l_n_v5)
            redu_Zagreb()
            #Second modified zagreb index
            def second_modi_Zagreb():
                sum6=0
                l6=[]
                l_n_v6=[]
                for q in range(1, h+1):
                    D6=f6[((h+1)*q)-(h+1)]
                    l6.append(D6)
                for i in l6:
                    sum6=sum6+i
                print("Second modified Zagreb index =",(sum6))
                for i in range(1,n_n_v+1):
                    n_v=sum6.subs(x,i) #n_v-numerical values
                    l_n_v6.append(n_v) 
                print('N.V of Second modified Zagreb index = ', l_n_v6)
            second_modi_Zagreb()
            #Randic index
            def Randic():
                sum7=0
                l7=[]
                l_n_v7=[]
                for q in range(1, h+1):
                    D7=f7[((h+1)*q)-(h+1)]
                    l7.append(D7)
                for i in l7:
                    sum7=sum7+i
                print("Randic index",(sum7))
                for i in range(1,n_n_v+1):
                    n_v=sum7.subs(x,i) #n_v-numerical values
                    l_n_v7.append(n_v) 
                print('N.V of Second modified Zagreb index = ', l_n_v7)
            Randic()
            #Reciprocal Randic index
            def reci_Randic():    
                sum8=0
                l8=[]
                l_n_v8=[]
                for q in range(1, h+1):
                    D8=f8[((h+1)*q)-(h+1)]
                    l8.append(D8)
                for i in l8:
                    sum8=sum8+i
                print("Reciprocal Randic index =",(sum8))
                for i in range(1,n_n_v+1):
                    n_v=sum8.subs(x,i) #n_v-numerical values
                    l_n_v8.append(n_v) 
                print('N.V of Reciprocal Randic index = ', l_n_v8)
            reci_Randic()
            #Reduced reciprocal Randic index
            def redu_reci_Randic():        
                sum9=0
                l9=[]
                l_n_v9=[]
                for q in range(1, h+1):
                    D9=f9[((h+1)*q)-(h+1)]
                    l9.append(D9)
                for i in l9:
                    sum9=sum9+i
                print("Reduced reciprocal Randic index =",(sum9))
                for i in range(1,n_n_v+1):
                    n_v=sum9.subs(x,i) #n_v-numerical values
                    l_n_v9.append(n_v) 
                print('N.V of Reduced reciprocal Randic index = ', l_n_v9)
            redu_reci_Randic()
            #General Randic index
            def gen_Randic():
                sum10=0
                l10=[]
                l_n_v10=[]
                for q in range(1, h+1):
                    D10=f10[((h+1)*q)-(h+1)]
                    l10.append(D10)
                for i in l10:
                    sum10=sum10+i
                print("General Randic index =",(sum10))
                for i in range(1,n_n_v+1):
                    n_v=sum10.subs(x,i) #n_v-numerical values
                    l_n_v10.append(n_v) 
                print('N.V of General Randic index = ', l_n_v10)
            gen_Randic()
            #Atom bond connectivity index
            def atm_bnd_connect():
                sum11=0
                l11=[]
                l_n_v11=[]
                for q in range(1, h+1):
                    D11=f11[((h+1)*q)-(h+1)]
                    l11.append(D11)
                for i in l11:
                    sum11=sum11+i
                print("Atom bond connectivity index =",(sum11))
                for i in range(1,n_n_v+1):
                    n_v=sum11.subs(x,i) #n_v-numerical values
                    l_n_v11.append(n_v) 
                print('N.V of Atom bond connectivity index = ', l_n_v11)
            atm_bnd_connect()
            #Geometric Arithematic index
            def geom_arith():
                sum12=0
                l12=[]
                l_n_v12=[]
                for q in range(1, h+1):
                    D12=f12[((h+1)*q)-(h+1)]
                    l12.append(D12)
                for i in l12:
                    sum12=sum12+i
                print("Geometric Arithematic index =",(sum12))
                for i in range(1,n_n_v+1):
                    n_v=sum12.subs(x,i) #n_v-numerical values
                    l_n_v12.append(n_v) 
                print('N.V of Geometric Arithematic index = ', l_n_v12)
            geom_arith()
            #Harmonic index
            def harm():            
                sum13=0
                l13=[]
                l_n_v13=[]
                for q in range(1, h+1):
                    D13=f13[((h+1)*q)-(h+1)]
                    l13.append(D13)
                for i in l13:
                    sum13=sum13+i
                print("Harmonic index =",(sum13))
                for i in range(1,n_n_v+1):
                    n_v=sum13.subs(x,i) #n_v-numerical values
                    l_n_v13.append(n_v) 
                print('N.V of Harmonic index = ', l_n_v13)
            harm()
            #Sum-connectivity index
            def sum_connect():
                sum14=0
                l14=[]
                l_n_v14=[]
                for q in range(1, h+1):
                    D14=f14[((h+1)*q)-(h+1)]
                    l14.append(D14)
                for i in l14:
                    sum14=sum14+i
                print("Sum-connectivity index =",(sum14))
                for i in range(1,n_n_v+1):
                    n_v=sum14.subs(x,i) #n_v-numerical values
                    l_n_v14.append(n_v) 
                print('N.V of Sum-connectivity index = ', l_n_v14)
            sum_connect()
            #Inverse sum index
            def inv_sum():
                sum15=0
                l15=[]
                l_n_v15=[]
                for q in range(1, h+1):
                    D15=f15[((h+1)*q)-(h+1)]
                    l15.append(D15)
                for i in l15:
                    sum15=sum15+i
                print("Inverse sum index =",(sum15))
                for i in range(1,n_n_v+1):
                    n_v=sum15.subs(x,i) #n_v-numerical values
                    l_n_v15.append(n_v) 
                print('N.V of Inverse sum index = ', l_n_v15)
            inv_sum()
            #Alberston index
            def alberst():
                sum16=0
                l16=[]
                l_n_v16=[]
                for q in range(1, h+1):
                    D16=f16[((h+1)*q)-(h+1)]
                    l16.append(D16)
                for i in l16:
                    sum16=sum16+i
                print("Alberston index =",(sum16))
                for i in range(1,n_n_v+1):
                    n_v=sum16.subs(x,i) #n_v-numerical values
                    l_n_v16.append(n_v) 
                print('N.V of Alberston index = ', l_n_v16)
            alberst()
            #Symmetric division index
            def symm_div():       
                sum17=0
                l17=[]
                l_n_v17=[]
                for q in range(1, h+1):
                    D17=f17[((h+1)*q)-(h+1)]
                    l17.append(D17)
                for i in l17:
                    sum17=sum17+i
                print("Symmetric division index =",(sum17))
                for i in range(1,n_n_v+1):
                    n_v=sum17.subs(x,i) #n_v-numerical values
                    l_n_v17.append(n_v) 
                print('N.V of Symmetric division index = ', l_n_v17)
            symm_div()
            #Forgotten index
            def forgot():
                sum18=0
                l18=[]
                l_n_v18=[]
                for q in range(1, h+1):
                    D18=f18[((h+1)*q)-(h+1)]
                    l18.append(D18)
                for i in l18:
                    sum18=sum18+i
                print("Forgotten index =",(sum18))
                for i in range(1,n_n_v+1):
                    n_v=sum18.subs(x,i) #n_v-numerical values
                    l_n_v18.append(n_v) 
                print('N.V of Forgotten index = ', l_n_v18)
            forgot() 
            #Sombor index
            def somb():
                sum19=0
                l19=[]
                l_n_v19=[]
                for q in range(1, h+1):
                    D19=f19[((h+1)*q)-(h+1)]
                    l19.append(D19)
                for i in l19:
                    sum19=sum19+i
                print("Sombor index =",(sum19))
                for i in range(1,n_n_v+1):
                    n_v=sum19.subs(x,i) #n_v-numerical values
                    l_n_v19.append(n_v) 
                print('N.V of Sombor index = ', l_n_v19)
            somb()
            #bi-Zagreb index
            def biZagreb():
                sum20=0
                l20=[]
                l_n_v20=[]
                for q in range(1, h+1):
                    D20=f20[((h+1)*q)-(h+1)]
                    l20.append(D20)
                for i in l20:
                    sum20=sum20+i
                print("bi-Zagreb index =",(sum20))
                for i in range(1,n_n_v+1):
                    n_v=sum20.subs(x,i) #n_v-numerical values
                    l_n_v20.append(n_v) 
                print('N.V of bi-Zagreb index = ', l_n_v20)
            biZagreb()
            #Tri-Zagreb index
            def TriZagreb():
                sum21=0
                l21=[]
                l_n_v21=[]
                for q in range(1, h+1):
                    D21=f21[((h+1)*q)-(h+1)]
                    l21.append(D21)
                for i in l21:
                    sum21=sum21+i
                print("Tri-Zagreb index =",(sum21))
                for i in range(1,n_n_v+1):
                    n_v=sum21.subs(x,i) #n_v-numerical values
                    l_n_v21.append(n_v) 
                print('N.V of Tri-Zagreb index = ', l_n_v21)
            TriZagreb()
            #Geometric Harmonic index
            def GeomHarm():
                sum22=0
                l22=[]
                l_n_v22=[]
                for q in range(1, h+1):
                    D22=f22[((h+1)*q)-(h+1)]
                    l22.append(D22)
                for i in l22:
                    sum22=sum22+i
                print("Geometric Harmonic index =",(sum22))
                for i in range(1,n_n_v+1):
                    n_v=sum22.subs(x,i) #n_v-numerical values
                    l_n_v22.append(n_v) 
                print('N.V of Geometric Harmonic index = ', l_n_v22)
            GeomHarm()
            #Geometric Bi-Zagreb index
            def GeomBiZagreb():
                sum23=0
                l23=[]
                l_n_v23=[]
                for q in range(1, h+1):
                    D23=f23[((h+1)*q)-(h+1)]
                    l23.append(D23)
                for i in l23:
                    sum23=sum23+i
                print("Geometric Bi-Zagreb index =",(sum23))
                for i in range(1,n_n_v+1):
                    n_v=sum23.subs(x,i) #n_v-numerical values
                    l_n_v23.append(n_v) 
                print('N.V of Geometric Bi-Zagreb index = ', l_n_v23)
            GeomBiZagreb()
            #Geometric Tri-Zagreb index
            def GeomTriZagreb():
                sum24=0
                l24=[]
                l_n_v24=[]
                for q in range(1, h+1):
                    D24=f24[((h+1)*q)-(h+1)]
                    l24.append(D24)
                for i in l24:
                    sum24=sum24+i
                print("Geometric Tri-Zagreb index =",(sum24))
                for i in range(1,n_n_v+1):
                    n_v=sum24.subs(x,i) #n_v-numerical values
                    l_n_v24.append(n_v) 
                print('N.V of Geometric Tri-Zagreb index = ', l_n_v24)
            GeomTriZagreb()
            #Harmonic-Geometric index
            def HarmGeom():
                sum25=0
                l25=[]
                l_n_v25=[]
                for q in range(1, h+1):
                    D25=f25[((h+1)*q)-(h+1)]
                    l25.append(D25)
                for i in l25:
                    sum25=sum25+i
                print("Harmonic-Geometric index =",(sum25))
                for i in range(1,n_n_v+1):
                    n_v=sum25.subs(x,i) #n_v-numerical values
                    l_n_v25.append(n_v) 
                print('N.V of Harmonic-Geometric index = ', l_n_v25)
            HarmGeom()
        except Exception:
            print("Invalid input detected")
    topoexpressd()
