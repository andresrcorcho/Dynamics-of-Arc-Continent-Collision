#!/usr/bin/env python
# coding: utf-8

# # 2D subduction-LayOut

# In[13]:


import UWGeodynamics as GEO
import numpy as np
import glucifer


# In[14]:


#Geometry functions- This is intented to be a fuctions module to evealuate models with distint geometries
import UWGeodynamics as GEO
import numpy as np


#Subducting plate creator
def SubductionCreator2D(Model,y0,thickness, dipAngle, dipLength, maxLength, orientation):
    #x and y are the (x0.,y0) coordinates where the shallow vertex of the slab will de constructed
    
    #Orientation= (-1)_ East dipping, (1)_ West dipping
    
    #Offset from model boundaries
    offset=100* u.kilometer
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
    #Trigonometry- Absolute dimensions of the dipping part
    angle = dipAngle * u.degree
    dx1 = np.cos(angle) * dipLength 
    dy1 = np.sin(angle) * dipLength 
    dx2 = np.sin(angle) * thickness 
    dy2 = np.cos(angle) * thickness 
    
    #Subducting plate calculation
    
    #Flat Segment Vertices    
    F1=(x0,y0)
    F2=(x0-(orientation*(maxLength-dipLength)),y0)
    F3=(x0-(orientation*(maxLength-dipLength)),y0-thickness)
    F4=(x0,y0-thickness)
    subducting2 = GEO.shapes.Polygon([F1,F2,F3,F4])
    #Dip Segment Vertices
    D1=(x0-(orientation*(maxLength-dipLength)),y0)
    D2=((x0-(orientation*(maxLength-dipLength)))-(orientation*dx1),(y0-dy1))
    D3=((x0-(orientation*(maxLength-dipLength)))-(orientation*(dx1-dx2)),(y0-dy1-dy2))
    D4=((x0-(orientation*(maxLength-dipLength)))+(orientation*dx2),(y0-dy2))
    subducting1 = GEO.shapes.Polygon([D1,D2,D3,D4])
    #Subducting plate
    subducting = subducting1 + subducting2
    #Calculation of Overriding plate geometry based in subducting plate geometry
    decoup=20* u.kilometer
    backd=400* u.kilometer
    #Flat Segment Vertices 
    f1=(X0,y0)
    f2=(((x0-(orientation*(maxLength-dipLength)))-(orientation*(thickness*dx1/dy1))),y0)
    f3=(((x0-(orientation*(maxLength-dipLength)))-(orientation*(thickness*dx1/dy1))-(orientation*decoup)),y0-thickness)
    f4=(X0,y0-thickness)
    over1= GEO.shapes.Polygon([f1,f2,f3,f4])
    
    #Dip Segment Vertices
    d1=(((x0-(orientation*(maxLength-dipLength)))-(orientation*(thickness*dx1/dy1))),y0)
    d2=(x0-(orientation*(maxLength-dipLength))-(orientation*decoup),y0)
    d3=(((x0-(orientation*(maxLength-dipLength)))-(orientation*(thickness*dx1/dy1))-(orientation*decoup)),y0-thickness)
    over2=GEO.shapes.Polygon([d1,d2,d3])
    #Overriding plate
    overriding = over1+over2
    
    #Decoupling air layer
    de1=(x0-(orientation*(maxLength-dipLength))-(orientation*decoup),y0)
    de2=(x0-(orientation*(maxLength-dipLength)),y0)
    de3=(((x0-(orientation*(maxLength-dipLength)))-(orientation*(thickness*dx1/dy1))),y0-thickness)
    de4=(((x0-(orientation*(maxLength-dipLength)))-(orientation*(thickness*dx1/dy1))-(orientation*decoup)),y0-thickness)
    
    decoup=GEO.shapes.Polygon([de1,de2,de3,de4])
    
    #Back Stop for oceanic plate
    b1=(x0,y0)
    b2=(x0-(orientation*backd),y0)
    b3=(x0-(orientation*backd),y0-(thickness/2))
    b4=(x0,y0-(thickness/2))
    
    backstop=GEO.shapes.Polygon([b1,b2,b3,b4])
    
    return subducting, overriding, decoup, backstop

#Overriding plate creator


# In[15]:


#Units
u = GEO.UnitRegistry


# In[16]:


#Scalling 1
#half_rate = 1.0 * u.centimeter / u.year
#model_length = 6000. * u.kilometer
#surfaceTemp = 273.15 * u.degK
#baseModelTemp = 1603.15 * u.degK
#bodyforce = 3300 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

#KL = model_length
#Kt = KL / half_rate
#KM = bodyforce * KL**2 * Kt**2
#KT = (baseModelTemp - surfaceTemp)

#GEO.scaling_coefficients["[length]"] = KL
#GEO.scaling_coefficients["[time]"] = Kt
#GEO.scaling_coefficients["[mass]"]= KM
#GEO.scaling_coefficients["[temperature]"] = KT


# In[21]:


#Scalling 2-- Stokes sink velocity
model_length = 1000. * u.kilometer
refDensity = 3300. * u.kilogram / u.meter**3
refViscosity = 1.4e19 * u.pascal * u.second

KL = model_length
KM = refDensity * KL**3
Kt = 1.0 / (refViscosity / KM * KL)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM


# In[22]:


#Model Dimensions
Model = GEO.Model(elementRes=(90, 90), 
                  minCoord=(-1500. * u.kilometer, -700.0 * u.kilometer), 
                  maxCoord=(1500. * u.kilometer, 0.0 * u.kilometer),
                  periodic = [False, False],
                  gravity=(0.0, -9.81 * u.meter / u.second**2))

#Model Output Folder
Model.outputDir = "outputs_SubductionOne"

#Model Thermodynamics
#Model.diffusivity = 1e-6 * u.metre**2 / u.second 
#Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)


# In[23]:


#Calculate Geometry
geometry=SubductionCreator2D(Model,0,100,30,300,1700,1)
#Materials in the model
#stickyAir = Model.add_material(name="Air", shape=GEO.shapes.Layer(top=Model.top, bottom=0. * u.kilometer))
Mantle = Model.add_material(name="Mantle", shape=GEO.shapes.Layer(top=0.*u.kilometer, bottom=Model.bottom))
OLithosphere = Model.add_material(name="Subducting Plate", shape=geometry[0])
Clithosphere= Model.add_material(name="Overriding plate", shape=geometry[1])
#decoup= Model.add_material(name="Decoupling", shape=geometry[2])
#backstop=Model.add_material(name="Back Stop", shape=geometry[3])


# In[24]:


#Preview of 2D materials-Materials Field (from swarm)
Fig = glucifer.Figure(figsize=(1200,400))
Fig.Points(Model.swarm, Model.materialField,fn_size=2.0, discrete=True,colours='blue white orange green white green')
Fig.show()


# In[10]:


#Tracers


# In[11]:


#Tracers Plot


# In[ ]:





# In[12]:


#Viscosities Registry
rh = GEO.ViscousCreepRegistry()
#Viscosity Setting
#stickyAir.viscosity = 1e19 * u.pascal * u.second
#decoup.viscosity= 1.4e19 * u.pascal * u.second
Mantle.viscosity =  1.4e19 * u.pascal * u.second#rh.Karato_and_Wu_1990# #
OLithosphere.visco sity = 200*1.4e21 * u.pascal * u.second
#Clithosphere.viscosity= 1.4e19 * u.pascal * u.second#1.4e23 * u.pascal * u.second#rh.Wang_et_al_2012 #
#backstop.viscosity=1e30 * u.pascal * u.second


# In[ ]:


#Plasticity Registry
#pl = GEO.PlasticityRegistry()
#Plasticity Setting
#OLithosphere.plasticity = pl.Huismans_et_al_2011_Crust
Clithosphere.plasticity = pl.Huismans_et_al_2011_Crust
#Mantle.plasticity = pl.Huismans_et_al_2011_Crust
OLithosphere.plasticity = GEO.VonMises(cohesion=48. * u.megapascal)
#Clithosphere.plasticity = GEO.VonMises(cohesion=48. * u.megapascal)
#decoup.plasticity=GEO.VonMises(cohesion=5. * u.megapascal)


# In[ ]:


#Density
#stickyAir.density = 1. * u.kilogram / u.metre**3
#decoup.density=3300. * u.kilogram / u.metre**3
Mantle.density = 3300. * u.kilogram / u.metre**3
OLithosphere.density =3400. * u.kilogram / u.metre**3
Clithosphere.density=2700. * u.kilogram / u.metre**3#2600. * u.kilogram / u.metre**3
#backstop.density=2700. * u.kilogram / u.metre**3


# In[ ]:


#Density Field
Fig = glucifer.Figure(figsize=(1200,400))
Fig.Points(Model.swarm, GEO.Dimensionalize(Model.densityField, u.kilogram / u.metre**3))
Fig.show()


# In[ ]:


#Viscosity Field
Fig = glucifer.Figure(figsize=(1200,400))
Fig.Points(Model.swarm, GEO.Dimensionalize(Model.viscosityField, u.pascal * u.second))
Fig.show()


# In[24]:


#Boundary Velocities- Subduction here is automatic.
#Model.set_velocityBCs(left=[0., None],
#                      right=[-5.0 * u.centimeter / u.year, None],
#                      bottom=[None, 0.],
#                      top=[None, 0.])

Model.set_velocityBCs(left=[0., None],
                      right=[0.,None],
                      bottom=[0., 0.],
                      top=[None, 0.])


# In[25]:


#Model Inicialization
Model.init_model()


# In[26]:


#Solver Parameters
GEO.rcParams["solver"] = "mumps"
GEO.rcParams["penalty"] = 1e6
GEO.rcParams["initial.nonlinear.tolerance"] = 1e-4


# In[ ]:


#Running Model
#!rm -rf output_FirstSubduction/
Model.run_for(duration=40*u.megayear, dt=0.5*u.megayear, checkpoint_interval=0.1*u.megayear,restartStep=False)
#Model.run_for(nstep=10, checkpoint_interval=0.1*u.megayear,restartStep=False)


# In[ ]:


#Materials at the end of the model
Fig = glucifer.Figure(figsize=(1200,400))
Fig.Points(Model.swarm, Model.materialField)
Fig.show()


# In[ ]:


#VIscosity Field at the end of the model
Fig = glucifer.Figure(figsize=(1200,400), title="Viscosity Field (Pa.s)", quality=3)
Fig.Points(Model.swarm, 
           GEO.Dimensionalize(Model.viscosityField, u.pascal * u.second),
           logScale=True,
           fn_size=3.0)
Fig.VectorArrows(Model.mesh, Model.velocityField)
Fig.show()


# In[ ]:


Fig = glucifer.Figure(figsize=(1200,400), title="Strain Rate Field (s-1)", quality=3)
Fig.Points(Model.swarm, 
           GEO.Dimensionalize(Model.strainRateField, 1.0/u.second),
           logScale=True,
           fn_size=3.0)
Fig.VectorArrows(Model.mesh, Model.velocityField)
Fig.show()


# In[ ]:


Fig = glucifer.Figure(figsize=(1200,400), title="Density Field (kg/cm^3)", quality=3)
Fig.Points(Model.swarm, 
           GEO.Dimensionalize(Model.densityField, u.kilogram/ u.cubic_centimeter),
           logScale=True,
           fn_size=3.0)
Fig.VectorArrows(Model.mesh, Model.velocityField)
Fig.show()


# # Parallel Run

# In[1]:


#Convert Script to .py

get_ipython().system('jupyter nbconvert --to=python MyFirstModel.ipynb')


# In[ ]:


#Run Script

get_ipython().system('mpirun -np 4 python MyFirstModel.py')


# In[ ]:




