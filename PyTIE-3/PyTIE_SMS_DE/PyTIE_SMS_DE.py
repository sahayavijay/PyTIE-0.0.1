class sms_topological_indices_Degree:
    def smstopoexpressd():
        try:
            print("Provide your data to calculate degree-based topological descriptors for the structure of a single molecule:")
            import math
            a1=[]
            a2=[]
            a3=[]
            a4=[]
            a5=[]
            a6=[]
            a7=[]
            a8=[]
            a9=[]
            a10=[]
            a11=[]
            a12=[]
            a13=[]
            a14=[]
            a15=[]
            a16=[]
            a17=[]
            a18=[]
            a19=[]
            a20=[]
            a21=[]
            a22=[]
            a23=[]
            a24=[]
            a25=[]
            b1=[]
            b2=[]
            b3=[]
            b6=[]
            b7=[]
            b8=[]
            b9=[]
            b10=[]
            b11=[]
            b12=[]
            b13=[]
            b14=[]
            b15=[]
            b16=[]
            b17=[]
            b18=[]
            b19=[]
            b20=[]
            b21=[]
            b22=[]
            b23=[]
            b24=[]
            b25=[]
            #ranges of edge partition
            m=int(input("Specify the minimum degree of a given graph (G): ")) 
            n=int(input("Specify the maximum degree of a given graph (G): "))
            for i in range(m, n+1):
                for j in range(m, n+1): 
                    if(i>=j):
                        k=int(input("Specify the values for partitioning edges in the graph of"+str((j, i))+":"))
                        #First Zagreb index
                        x1=k*(j+i)
                        a1.append(x1)
                        #First Zagreb entropy
                        y1=k*((j+i)*math.log(j+i))
                        b1.append(y1)
                        #Second Zagreb index
                        x2=k*(j*i)
                        a2.append(x2)
                        #Second Zagreb entropy
                        y2=k*((j*i)*math.log(j*i))
                        b2.append(y2)
                        #Hyper Zagreb index
                        x3=k*((j+i)**2)
                        a3.append(x3)
                        #Hyper Zagreb entropy
                        y3=k*(((j+i)**2)*math.log((j+i)**2))
                        b3.append(y3)
                        #Third Zagreb index
                        x4=k*(abs(j-i)) 
                        a4.append(x4)
                        #Reduced Zagreb index
                        x5=k*((j-1)*(i-1))
                        a5.append(x5)
                        #Second modified Zagreb index
                        x6=k*(1/(j*i)) 
                        a6.append(x6)
                        #Second modified Zagreb entropy
                        y6=k*((1/(j*i))*math.log((1/(j*i))))
                        b6.append(y6)
                        #Randic index
                        x7=k*(1/math.sqrt(j*i)) 
                        a7.append(x7)
                        #Randic entropy
                        y7=k*((1/math.sqrt(j*i))*math.log((1/math.sqrt(j*i))))
                        b7.append(y7)
                        #Reciprocal Randic index
                        x8=k*(math.sqrt(j*i))
                        a8.append(x8)
                        #Reciprocal Randic entropy
                        y8=k*((math.sqrt(j*i))*math.log((math.sqrt(j*i))))
                        b8.append(y8)
                        #Reduced reciprocal Randic index
                        x9=k*(math.sqrt((j-1)*(i-1)))
                        a9.append(x9)
                        #General Randic index
                        x10=k*((j*i)**(-1)) 
                        a10.append(x10)
                        #General Randic entropy
                        y10=k*(((j*i)**(-1))*math.log(((j*i)**(-1))))
                        b10.append(y10)
                        #Atom bond connectivity index
                        x11=k*(math.sqrt((j+i-2)/(j*i))) 
                        a11.append(x11)
                        #Geometric Arithematic index
                        x12=k*(2*(math.sqrt(j*i)/(j+i)))
                        a12.append(x12)
                        #Geometric Arithematic entropy
                        y12=k*((2*(math.sqrt(j*i)/(j+i)))*math.log((2*(math.sqrt(j*i)/(j+i)))))
                        b12.append(y12)
                        #Harmonic index
                        x13=k*(2/(j+i))
                        a13.append(x13)
                        #Harmonic entropy
                        y13=k*((2/(j+i))*math.log((2/(j+i))))
                        b13.append(y13)
                        #Sum-connectivity index
                        x14=k*(1/math.sqrt(j+i))
                        a14.append(x14)
                        #Sum-connectivity entropy
                        y14=k*((1/math.sqrt(j+i))*math.log((1/math.sqrt(j+i))))
                        b14.append(y14)
                        #Inverse sum index
                        x15=k*((j*i)/(j+i))
                        a15.append(x15)
                        #Inverse sum entropy
                        y15=k*(((j*i)/(j+i))*math.log(((j*i)/(j+i))))
                        b15.append(y15)
                        #Alberston index
                        x16=k*(abs(j-i)) 
                        a16.append(x16)
                        #Symmetric division index
                        x17=k*(((j**2)+(i**2))/(j*i)) 
                        a17.append(x17)
                        #Symmetric division entropy
                        y17=k*((((j**2)+(i**2))/(j*i))*math.log((((j**2)+(i**2))/(j*i))))
                        b17.append(y17)
                        #Forgotten index
                        x18=k*((j**2)+(i**2))
                        a18.append(x18)
                        #Forgotten entropy
                        y18=k*(((j**2)+(i**2))*math.log(((j**2)+(i**2))))
                        b18.append(y18)
                        #Sombor index
                        x19=k*(math.sqrt((j**2)+(i**2)))
                        a19.append(x19)
                        #Sombor entropy
                        y19=k*((math.sqrt((j**2)+(i**2)))*math.log((math.sqrt((j**2)+(i**2)))))
                        b19.append(y19)
                        #bi-Zagreb index
                        x20=k*(j+i+(j*i))
                        a20.append(x20)
                        #bi-Zagreb entropy
                        y20=k*((j+i+(j*i))*math.log((j+i+(j*i))))
                        b20.append(y20)
                        #Tri-Zagreb index
                        x21=k*((j**2)+(i**2)+(j*i))
                        a21.append(x21)
                        #Tri-Zagreb entropy
                        y21=k*(((j**2)+(i**2)+(j*i))*math.log(((j**2)+(i**2)+(j*i))))
                        b21.append(y21)
                        #Geometric Harmonic index
                        x22=k*((math.sqrt(j*i)*(j+i))/2)
                        a22.append(x22)
                        #Geometric Harmonic entropy
                        y22=k*((((math.sqrt(j*i)*(j+i))/2)*math.log((math.sqrt(j*i)*(j+i))/2)))
                        b22.append(y22)
                        #Geometric Bi-Zagreb index
                        x23=k*(math.sqrt(j*i)/(j+i+j*i))
                        a23.append(x23)
                        #Geometric Bi-Zagreb entropy
                        y23=k*((math.sqrt(j*i)/(j+i+j*i))*math.log((math.sqrt(j*i)/(j+i+j*i))))
                        b23.append(y23)
                        #Geometric Tri-Zagreb index
                        x24=k*(math.sqrt(j*i)/((j**2)+(i**2)+j*i))
                        a24.append(x24)
                        #Geometric Tri-Zagreb entropy
                        y24=k*((math.sqrt(j*i)/((j**2)+(i**2)+j*i))*math.log((math.sqrt(j*i)/((j**2)+(i**2)+j*i))))
                        b24.append(y24)
                        #Harmonic Geometric index
                        x25=k*(2/((math.sqrt(j*i))*(j+i)))
                        a25.append(x25)
                        #Harmonic Geometric entropy
                        y25=k*((2/((math.sqrt(j*i))*(j+i)))*math.log((2/((math.sqrt(j*i))*(j+i)))))
                        b25.append(y25)
            sum=0
            for i in a1:
                c1=sum=sum+i
            print("First Zagreb index=",c1)
            sum=0
            for i in a2:
                c2=sum=sum+i
            print("Second Zagreb index=",c2)
            sum=0
            for i in a3:
                c3=sum=sum+i
            print("Hyper Zagreb index=",c3)
            sum=0
            for i in a4:
                c4=sum=sum+i
            print("Third Zagreb index=",c4)
            sum=0
            for i in a5:
                c5=sum=sum+i
            print("Reduced Zagreb index=",c5)
            sum=0
            for i in a6:
                c6=sum=sum+i
            print("Second modified Zagreb index=",sum)
            sum=0
            for i in a7:
                c7=sum=sum+i
            print("Randic index=",c7)
            sum=0
            for i in a8:
                c8=sum=sum+i
            print("Reciprocal Randic index=",c8)
            sum=0
            for i in a9:
                c9=sum=sum+i
            print("Reduced reciprocal Randic index=",c9)
            sum=0
            for i in a10:
                c10=sum=sum+i
            print("General Randic index=",c10)
            sum=0
            for i in a11:
                c11=sum=sum+i
            print("Atom bond connectivity index=",c11)
            sum=0
            for i in a12:
                c12=sum=sum+i
            print("Geometric Arithematic index=",c12)
            sum=0
            for i in a13:
                c13=sum=sum+i
            print("Harmonic index=",c13)
            sum=0
            for i in a14:
                c14=sum=sum+i
            print("Sum-connectivity index=",c14)
            sum=0
            for i in a15:
                c15=sum=sum+i
            print("Inverse sum index=",c15)
            sum=0
            for i in a16:
                c16=sum=sum+i
            print("Alberston index=",c16)
            sum=0
            for i in a17:
                c17=sum=sum+i
            print("Symmetric division index=",c17)
            sum=0
            for i in a18:
                c18=sum=sum+i
            print("Forgotten index=",c18)
            sum=0
            for i in a19:
                c19=sum=sum+i
            print("Sombor index=",c19)
            sum=0
            for i in a20:
                c20=sum=sum+i
            print("bi-Zagreb index=",c20)
            sum=0
            for i in a21:
                c21=sum=sum+i
            print("Tri-Zagreb index=",c21)
            sum=0
            for i in a22:
                c22=sum=sum+i
            print("Geometric Harmonic index =",c22)
            sum=0
            for i in a23:
                c23=sum=sum+i
            print("Geometric Bi-Zagreb index =",c23)
            sum=0
            for i in a24:
                c24=sum=sum+i
            print("Geometric Tri-Zagreb index =",c24)
            sum=0
            for i in a25:
                c25=sum=sum+i
            print("Harmonic Geometric index =",c25)
            sum=0
            for i in b1:
                d1=sum=sum+i
            z1=math.log(c1)
            print("First Zagreb entropy=",z1-((1/c1)*d1))
            sum=0
            for i in b2:
                d2=sum=sum+i
            z2=math.log(c2)
            print("Second Zagreb entropy=",z2-((1/c2)*d2))
            sum=0
            for i in b3:
                d3=sum=sum+i
            z3=math.log(c3)
            print("Hyper Zagreb entropy=",z3-((1/c3)*d3))
            sum=0
            for i in b6:
                d6=sum=sum+i
            z6=math.log(c6)
            print("Second modified Zagreb entropy=",z6-((1/c6)*d6))
            sum=0
            for i in b7:
                d7=sum=sum+i
            z7=math.log(c7)
            print("Randic entropy=",z7-((1/c7)*d7))
            sum=0
            for i in b8:
                d8=sum=sum+i
            z8=math.log(c8)
            print("Reciprocal Randic entropy=",z8-((1/c8)*d8))
            sum=0
            for i in b10:
                d10=sum=sum+i
            z10=math.log(c10)
            print("General Randic entropy=",z10-((1/c10)*d10))
            sum=0
            for i in b12:
                d12=sum=sum+i
            z12=math.log(c12)
            print("Geometric Arithematic entropy=",z12-((1/c12)*d12))
            sum=0
            for i in b13:
                d13=sum=sum+i
            z13=math.log(c13)
            print("Harmonic entropy=",z13-((1/c13)*d13))
            sum=0
            for i in b14:
                d14=sum=sum+i
            z14=math.log(c14)
            print("Sum-connectivity entropy=",z14-((1/c14)*d14))
            sum=0
            for i in b15:
                d15=sum=sum+i
            z15=math.log(c15)
            print("Inverse sum entropy=",z15-((1/c15)*d15))
            sum=0
            for i in b17:
                d17=sum=sum+i
            z17=math.log(c17)
            print("Symmetric division entropy=",z17-((1/c17)*d17))
            sum=0
            for i in b18:
                d18=sum=sum+i
            z18=math.log(c18)
            print("Forgotten entropy=",z18-((1/c18)*d18))
            sum=0
            for i in b19:
                d19=sum=sum+i
            z19=math.log(c19)
            print("Sombor entropy=",z19-((1/c19)*d19))
            sum=0
            for i in b20:
                d20=sum=sum+i
            z20=math.log(c20)
            print("bi-Zagreb entropy=",z20-((1/c20)*d20))
            sum=0
            for i in b21:
                d21=sum=sum+i
            z21=math.log(c21)
            print("Tri-Zagreb entropy=",z21-((1/c21)*d21))
            sum=0
            for i in b22:
                d22=sum=sum+i
            z22=math.log(c22)
            print("Geometric Harmonic entropy=",z22-((1/c22)*d22))
            sum=0
            for i in b23:
                d23=sum=sum+i
            z23=math.log(c23)
            print("Geometric Bi-Zagreb entropy=",z23-((1/c23)*d23))
            sum=0
            for i in b24:
                d24=sum=sum+i
            z24=math.log(c24)
            print("Geometric Geometric Tri-Zagreb entropy=",z24-((1/c24)*d24))
            sum=0
            for i in b25:
                d25=sum=sum+i
            z25=math.log(c25)
            print("Geometric Harmonic Geometric entropy=",z25-((1/c25)*d25))
        except Exception:
            print("Invalid input detected")
    smstopoexpressd()
sms_topological_indices_Degree()
