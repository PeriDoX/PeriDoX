import sys

#xMax = 16
#xMax = int(input("Anzahl in x: "))
#yMax = int(input("Anzahl in y: "))
#zMax = int(input("Anzahl in z: "))

L = 0.050
h = 0.020
B = 0.006
a = 0.005

xMax = 51
yMax = 22
zMax = 7

#dx = float(input("Schrittweite dx: "))
#dy = float(input("Schrittweite dy: "))
#dz = float(input("Schrittweite dz: "))

dx = L/(xMax-1)
dy = h/(yMax-1)
dz = B/(zMax-1)

print(dx,dy,dz)


delta = 1 # int(input("Randbereich: "))
i = 0
d = dx * dy * dz

#xMin = -(xMax/2 - 0.5) * dx
xMin = -a
yMin = -(yMax/2 - 0.5) * dy
zMin = -(zMax/2 - 0.5) * dz

plateFile = open("plate.txt","w")
nodeset_1 = open("nodeset_1.txt","w")
nodeset_2 = open("nodeset_2.txt","w")
nodeset_3 = open("nodeset_3.txt","w")
nodeset_4 = open("nodeset_4.txt","w")

plateFile.write("# x y z block_id volume\n")

for x in range(1, xMax+1):
    for y in range(1, yMax+1):
        for z in range(1, zMax+1):
            #print "%f %f %f" % (xMin + x, yMin + y, zMin + z)
            i += 1
            k = 1


            #if xMin + (x-1) * dx<-a/2:
            #  k=2;
            
            plateFile.write("%f %f %f %d %1.11f\n" % (xMin + (x-1) * dx, yMin + (y-1) * dy, zMin + (z-1) * dz, k, d))
            # Fig 1 one half of support span    
            if x < 4:
                if yMin + (y-1) * dy>0:
                
                    nodeset_1.write("%d\n" % (i))
            if x < 4:
                if yMin + (y-1) * dy<0:
                    nodeset_2.write("%d\n" % (i))
            if z < 2:
                
                    nodeset_3.write("%d\n" % (i))
            if z > zMax:
                
                    nodeset_3.write("%d\n" % (i))


plateFile.close()
nodeset_1.close()
nodeset_2.close()
nodeset_3.close()
nodeset_4.close()
