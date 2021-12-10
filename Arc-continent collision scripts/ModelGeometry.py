import UWGeodynamics as GEO
import numpy as np
from underworld import function as fn
#Units
u = GEO.UnitRegistry

#Subducting plate LayerCreator
def interpolateTracer(coord1,coord2,nPoints):
    x=np.linspace(coord1[0],coord2[0],nPoints)*u.kilometer
    m=(coord2[1]-coord1[1])/(coord2[0]-coord1[0])
    y=m*(x-coord1[0])+coord1[1]
    return x,y

def fuseInList(List):
    return List

def limitArray(ListI,maxAmount):
    ListF=[]
    counter=0
    for i in ListI:
        if counter<maxAmount:
            ListF.append(i)
        counter=counter+1
    return ListF
#Fuse list at the same dimensional level- But same Length
def fuseList(lists):
    fListx=[]
    fListy=[]
    for i in lists: 
        fauxX=[i[0]]
        fauxY=[i[1]]
        fListx=fListx+fauxX
        fListy=fListy+fauxY
    return fListx,fListy
#Multidimentional case
def fuseListM(lists):
    fListx=[]
    fListy=[]
    for i in lists: 
        for j,k in zip(i[0],i[1]):
#         fauxX=[i[0]]
#         fauxY=[i[1]]
            fListx=fListx+[j]
            fListy=fListy+[k]
    return fListx,fListy

def isRepeated(numTuple,List):
    ans=False
    counter=0
    ListX=List[0]
    ListY=List[1]
    for i,j in zip(ListX,ListY):
        if numTuple[0]==i and numTuple[1]==j:
            counter=counter+1
    if counter>1:
        ans=True
    return ans, counter

def isInList(numTuple,List):
    ans=False
    counter=0
    ListX=List[0]
    ListY=List[1]
    for i,j in zip(ListX,ListY):
        if numTuple[0]==i and numTuple[1]==j:
            counter=counter+1
            ans=True
    return ans,counter

def rmRepeated(ListI):
    listFx=[]
    listFy=[]
    ListX=ListI[0]
    ListY=ListI[1]
    counter=0
    for i,j in zip(ListX,ListY):
        if counter==0:
            listFx.append(i)
            listFy.append(j)
        else:
            if isInList((i,j),[listFx,listFy])[0]==False:
                listFx.append(i)
                listFy.append(j)
        counter=counter+1
    return listFx,listFy

def ListToNd(listI):
    ListF=[]
    for i in listI:
        ListF.append(GEO.nd(i))
    return ListF

def SubLayerCreator(x0,y0,Lthickness,dipAngle, dipLength, maxLength,endOff0,endOff, orientation,ExLen,BarcT,bStrip=False):
    #Trigonometry- Absolute dimensions of the dipping part
    angle = dipAngle * u.degree
    dx1 = np.cos(angle) * dipLength 
    dy1 = np.sin(angle) * dipLength 
    dx2 = np.sin(angle) * Lthickness 
    dy2 = np.cos(angle) * Lthickness 
    #Flat Segment Vertices    
    F1=(x0-(orientation*endOff0),y0)
    F2=(x0-(orientation*(maxLength-dipLength)),y0)
    F3=(x0-(orientation*(maxLength-dipLength)),y0-Lthickness)
    F4=(x0-(orientation*(endOff)),y0-Lthickness)
    subducting2 = GEO.shapes.Polygon([F1,F2,F3,F4])
    #Flat segment in sections
    ns=4
    lenRatio=maxLength/4
    conPLen=lenRatio-dipLength
    xCero=(x0-(orientation*(maxLength-dipLength))+ (conPLen*orientation))
    seg=[]
    for i in range(0,ns-2):
        vx1=(xCero,y0)
        vx2=((xCero+(orientation*lenRatio)),y0)
        vx3=((xCero+(orientation*lenRatio)),y0-Lthickness)
        vx4=(xCero,y0-Lthickness)
        seg.append(GEO.shapes.Polygon([vx1,vx2,vx3,vx4]))
        xCero=xCero+(orientation*lenRatio)
    
    #Dip Segment Vertices
    D1=(x0-(orientation*(maxLength-dipLength)),y0)
    D2=((x0-(orientation*(maxLength-dipLength)))-(orientation*dx1),(y0-dy1))
    D3=((x0-(orientation*(maxLength-dipLength)))-(orientation*dx1),(y0-dy1-Lthickness))
    D4=((x0-(orientation*(maxLength-dipLength)))-(orientation*(dx1-dx2)),(y0-dy1-dy2))
    D5=((x0-(orientation*(maxLength-dipLength)))+(orientation*dx2),(y0-dy2))
    subducting1 = GEO.shapes.Polygon([D1,D2,D3,D4,D5])
    subducting = subducting1 | subducting2
    SlabTracers=fuseList([F1,D1,D2,D3,F3,F4])
    #As this function only return a layer, it could include a layer of exotic material/strip with vayring lenght, 
    #which will be located always in the same location of the layer (out of the dipping segment)
    #Subducting plate
    #Exotic crust Segment
    exD=450* u.kilometer #Distance of Exotic terrane
    Exlen=ExLen
    ExLen=ExLen* u.kilometer
    wlen=BarcT* u.kilometer
    arcL=0
    #Weak layer
    we1=((x0-(orientation*(maxLength-dipLength)))+(orientation*(exD+ExLen)),y0)
    we2=((x0-(orientation*(maxLength-dipLength)))+(orientation*(exD+ExLen+wlen)),y0)
    we3=((x0-(orientation*(maxLength-dipLength)))+(orientation*(exD+ExLen+wlen)),y0-Lthickness)
    we4=((x0-(orientation*(maxLength-dipLength)))+(orientation*(exD+ExLen)),y0-Lthickness)
    weakk=GEO.shapes.Polygon([we1,we2,we3,we4])
    
    if Exlen>0:
        sa1=((x0-(orientation*(maxLength-dipLength)))+(orientation*exD),y0)
        sa2=((x0-(orientation*(maxLength-dipLength)))+(orientation*(exD+ExLen)),y0)
        sa3=((x0-(orientation*(maxLength-dipLength)))+(orientation*(exD+ExLen)),y0-Lthickness)
        sa4=((x0-(orientation*(maxLength-dipLength)))+(orientation*exD),y0-Lthickness)
        arcA=GEO.shapes.Polygon([sa1,sa2,sa3,sa4])
        arcL=arcA   
    ArcTracers=fuseList([sa1,sa2,sa3,sa4])
    #Buoyant Strip
    bStripL=0
    buD=lenRatio
    if bStrip==True:
        bs1=(x0-(orientation*endOff0),y0)
        bs2=(x0-(orientation*buD),y0)
        bs3=(x0-(orientation*buD),y0-Lthickness)
        bs4=(x0-(orientation*endOff),y0-Lthickness)
        bStripL=GEO.shapes.Polygon([bs1,bs2,bs3,bs4])
        
    return [subducting1,subducting2], arcL, weakk,SlabTracers,ArcTracers,bStripL,seg

def OverCreatorL(X0,x0,y0,thickness,decoup,dipAngle, dipLength, maxLength, orientation,offs,Lprop,Tprop):
    #Calculation of Overriding plate geometry based in subducting plate geometry
    #Trigonometry- Absolute dimensions of the dipping part
    angle = dipAngle * u.degree
    dx1 = np.cos(angle) * dipLength 
    dy1 = np.sin(angle) * dipLength 
    dx2 = np.sin(angle) * thickness 
    dy2 = np.cos(angle) * thickness 
    #Flat Segment Vertices 
    f1=(X0,y0)
    f2=(((x0-(orientation*(maxLength-dipLength)))-(orientation*np.absolute(((x0-(orientation*(maxLength-dipLength)))-(X0))/Lprop))),y0)
    f3=(((x0-(orientation*(maxLength-dipLength)))-(orientation*np.absolute(((x0-(orientation*(maxLength-dipLength)))-(X0))/Lprop))),y0-thickness)
    f4=(X0,y0-thickness)
    over1= GEO.shapes.Polygon([f1,f2,f3,f4])
    TracersCraton=fuseList([f1,f2,f3,f4])
    #Overriding plate
    overriding = over1
    #Correction of y0- As weak OP is less thick
    if y0<0 * u.kilometer:
        y0=y0+(thickness-((thickness/Tprop)))
    #Weak OP
    p1=(((x0-(orientation*(maxLength-dipLength)))-(orientation*np.absolute(((x0-(orientation*(maxLength-dipLength)))-(X0))/Lprop))),y0)
    p2=(x0-(orientation*(maxLength-dipLength))-(orientation*decoup)-(orientation*offs),y0)
    p3=(((x0-(orientation*(maxLength-dipLength)))-(orientation*((thickness/Tprop)*dx1/dy1))-(orientation*decoup)-(orientation*offs)),y0-(thickness/Tprop))
    p4=(((x0-(orientation*(maxLength-dipLength)))-(orientation*np.absolute(((x0-(orientation*(maxLength-dipLength)))-(X0))/Lprop))),y0-(thickness/Tprop))
    weak=GEO.shapes.Polygon([p1,p2,p3,p4]) 
    TracersWeak=fuseList([p1,p2,p3,p4])
    #Decoupling  layer
    thickness=thickness/2.
    de1=(x0-(orientation*(maxLength-dipLength))-(orientation*decoup)-(orientation*offs),y0)
    de2=(x0-(orientation*(maxLength-dipLength))-(orientation*offs),y0)
    de3=(((x0-(orientation*(maxLength-dipLength)))-(orientation*((thickness/Tprop)*dx1/dy1))-(orientation*offs)),y0-(thickness/Tprop))
    de4=(((x0-(orientation*(maxLength-dipLength)))-(orientation*((thickness/Tprop)*dx1/dy1))-(orientation*decoup)-(orientation*offs)),y0-(thickness/Tprop))
    decou1=GEO.shapes.Polygon([de1,de2,de3,de4])
    
    y0=y0-((thickness/Tprop))
    offsetL=(((thickness/Tprop)*dx1/dy1))
    offs=(offs+offsetL)
    De1=(x0-(orientation*(maxLength-dipLength))-(orientation*decoup)-(orientation*offs),y0)
    De2=(x0-(orientation*(maxLength-dipLength))-(orientation*offs),y0)
    De3=(((x0-(orientation*(maxLength-dipLength)))-(orientation*((thickness/Tprop)*dx1/dy1))-(orientation*offs)),y0-(thickness/Tprop))
    De4=(((x0-(orientation*(maxLength-dipLength)))-(orientation*((thickness/Tprop)*dx1/dy1))-(orientation*decoup)-(orientation*offs)),y0-(thickness/Tprop))
    decou2=GEO.shapes.Polygon([De1,De2,De3,De4])
    return overriding, weak, decou1, decou2,TracersWeak,TracersCraton

import UWGeodynamics as GEO
import numpy as np
#Subducting plate creator
def SubductionCreator(Model,y0,thickness, dipAngle, dipLength, maxLength, orientation, SLayers,OLayers, ExLens,BarcTs,bStrips,decoup):
    #x and y are the (x0.,y0) coordinates where the shallow vertex of the slab will de constructed
    #Orientation= (-1)_ East dipping, (1)_ West dipping
    #Offset from model boundaries
    offset=100.* u.kilometer
    angle = dipAngle * u.degree
    dx1 = np.cos(angle) * dipLength 
    dy1 = np.sin(angle) * dipLength
    #Units for parameters
    if orientation==1:
        x0=Model.maxCoord[0]-offset
        X0=Model.minCoord[0]+offset
    else:
        x0=Model.minCoord[0]+offset
        X0=Model.maxCoord[0]-offset
        
    y0=y0 * u.kilometer
    thickness= thickness * u.kilometer
    dipLength= dipLength * u.kilometer
    maxLength= maxLength * u.kilometer
    
    #Data from subducting plate subduction 
    #Subducting plate calculation
    subducting=[]
    subductingD=[]
    subductingF=[]
    arc=[]
    bstop=[]
    sections=[]
    weaak=[]
    ratio=thickness/ SLayers
    lay=0
    #eRat=(200.* u.kilometer)/ SLayers
    eRat=0.* u.kilometer
    eOff0=0.* u.kilometer #300* u.kilometer
    start=False
     #[subducting1_Dip,subducting2_Flat], arcL, bStripL,seg,weakk
    
    for i in range (0,SLayers):
        eOff=eOff0+(eRat)#eOff=eOff0#
        subductingData=SubLayerCreator(x0,y0-lay,ratio,dipAngle, dipLength, maxLength,eOff0,eOff,orientation,ExLens[i],BarcTs[i],bStrips[i])
        subductingD.append(subductingData[0][0])
        subductingF.append(subductingData[0][1])
        subducting.append(subductingData[0][0] | subductingData[0][1])
        weaak.append(subductingData[2])
        if ExLens[i]>0:
            arc.append(subductingData[1])
        if bStrips[i]==True:
            bstop.append(subductingData[5])
        sections.append(subductingData[6])
        auxT=subductingData[3]
        auxAT=subductingData[4]
        if start==False:
            SPTracers=auxT
            ArcTracers=auxAT
            start=True
        else:
            SPTracers=fuseListM([SPTracers,auxT])
            ArcTracers=fuseListM([ArcTracers,auxAT])
        lay=lay+ratio
        eOff0=eOff0+eRat
        
    #Sum All Arc and backstop polygons?? but maybe arc not. They can include different properties
    #Exotic crust calculation
    #Calculation of Overriding plate geometry based in subducting plate geometry
    decoup=decoup * u.kilometer
    sP=[]
    wP=[]
    dec=[]
    OPTracers=[]
    ratio2=(thickness/(10./15.))/OLayers
    dlay=0
    offset0L=0* u.kilometer
    offsetL=(((thickness/OLayers)*dx1/dy1))
    start=False
    
    for j in range (0,OLayers):
        ovedata=OverCreatorL(X0,x0,y0-dlay,ratio2,decoup,dipAngle, dipLength, maxLength, orientation,offset0L,Lprop=(1.5),Tprop=(15./10.))
        sP.append(ovedata[0])
        wP.append(ovedata[1])
        dec.append(ovedata[2])
        dec.append(ovedata[3])
        auxT=ovedata[4]
        cratonT=ovedata[5]
        if start==False:
            OPTracers=auxT
            CratTracers=cratonT
            start=True
        else:
            OPTracers=fuseListM([OPTracers,auxT])
            CratTracers=fuseListM([CratTracers,cratonT])
        dlay=dlay+ratio2
        offset0L=offset0L+offsetL
    #X coords in plates interface    
    xlimit=(x0-(orientation*(maxLength-dipLength)))
    #crLimit=((x0-(orientation*(maxLength-dipLength)))-(orientation*np.absolute(((x0-(orientation*(maxLength-dipLength)))-(X0))/Lprop)))
    
    
    return [subductingD,subductingF,subducting],sP, wP,arc,dec, weaak,[orientation,xlimit],SPTracers,OPTracers,ArcTracers,CratTracers,bstop,sections

