
//_____________________________G4BoundaryPhaseContrast.cc_________________________

#include "G4BoundaryPhaseContrast.hh"



//_________________________________FUNCTIONS______________________________


//1st function : dot product 
//this function is optional

G4double G4BoundaryPhaseContrast::dot(G4double d[],G4double n[]){

    G4double dx=d[0]*n[0];
    G4double dy=d[1]*n[1];
    G4double dz=d[2]*n[2];

    G4double prod= dx+dy+dz;
    G4double cosang=prod;
    G4double ang_rad= acos(cosang);
    G4double ang=ang_rad*180/3.1415;
    return ang;
    }

//2nd function = Calculating the gradient (normal vector)

G4double G4BoundaryPhaseContrast::diff (G4double a[],G4double x,G4double y,G4double z, G4int n){   //..................................................// 
    G4double dif=0;										   //............. a=fitted parameters.................//
    if (n==1){											   //............. n=1 for x ,=2 for y ,=3 for z.......//
        dif=2*a[0]*x + a[3]*y+ a[4]*z+a[6];							   //..................................................//
    }
    if(n==2){
        dif= 2*a[1]*y + a[3]*x +a[5]*z+a[7];
    }
    if(n==3){
        dif=2*a[2]*z +a[4]*x+a[5]*y+a[8];
    }

    return dif;
    }


//3rd function = find boundary
//fitted ellipsoidal function 

G4double G4BoundaryPhaseContrast::boundary(G4double a[],G4double x,G4double y ,G4double z){

    G4double fitted = a[0]* pow(x,2)+ a[1]*pow(y,2)+ a[2]*pow(z,2)+ a[3]*x*y +a[4]*x*z +a[5]*y*z + a[6]*x + a[7]*y + a[8]*z +a[9];

    G4double marz=0;
    if (fitted < .0005 && fitted > -.0005){
        marz=1;  
    }
    else{
        marz=0;
    }
    return marz;
}


//4th function = which fitted function ?

G4double G4BoundaryPhaseContrast::fit(G4double x, G4double y , G4double z){			//.....................................................................................//
G4double f=0;											//......outputs of the python function (fitting function) should be copied here!.......//
 if (x>=0 && y>=0 && z>=0 ){f=1;}								//.....................................................................................//
 if (x>0 && y>0 && z<=0 ){f=2;}
 if (x>0 && y<0 && z>=0 ){f=3;}
 if (x<0 && y>0 && z>=0 ){f=4;}
 if (x<0 && y<0 && z>=0 ){f=5;}
 if (x<0 && y>0 && z<=0 ){f=6;}
 if (x>0 && y<0 && z<=0 ){f=7;}
 if (x<0 && y<0 && z<=0 ){f=8;}
return f;
}
void G4BoundaryPhaseContrast::index(G4double f,G4double a[]){
if (f==1){
a[0] = 62492.42674610019 ;
a[1] = 51037.84087750316 ;
a[2] = 822233.4155580997 ;
a[3] = -63218.53164168447 ;
a[4] = 32078.110543847084 ;
a[5] = -72749.28930270672 ;
a[6] = -30.01706966315396 ;
a[7] = 24.817077849293128 ;
a[8] = 213.92182709183544 ;
a[9] = -1.0 ;
}
if (f==2){
a[0] = 62606.78170695901 ;
a[1] = 51136.918440297246 ;
a[2] = 832667.3562214375 ;
a[3] = -63346.88473594934 ;
a[4] = -33744.77237933874 ;
a[5] = 72719.78481730819 ;
a[6] = -30.374276688322425 ;
a[7] = 24.862324570363853 ;
a[8] = -201.99640310788527 ;
a[9] = -1.0 ;
}
if (f==3){
a[0] = 62492.42674610019 ;
a[1] = 51037.84087750316 ;
a[2] = 822233.4155580997 ;
a[3] = 63218.53164168447 ;
a[4] = 32078.110543847084 ;
a[5] = 72749.28930270672 ;
a[6] = -30.01706966315396 ;
a[7] = -24.817077849293128 ;
a[8] = 213.92182709183544 ;
a[9] = -1.0 ;
}
if (f==4){
a[0] = 55097.694467663765 ;
a[1] = 45237.90702147782 ;
a[2] = 734885.987244606 ;
a[3] = 61271.89123310894 ;
a[4] = -24113.774471193552 ;
a[5] = -81251.17187865078 ;
a[6] = 2.6013991633662954 ;
a[7] = 47.47156475123484 ;
a[8] = 267.2414241330698 ;
a[9] = -1.0 ;
}
if (f==5){
a[0] = 55097.694467663765 ;
a[1] = 45237.90702147782 ;
a[2] = 734885.987244606 ;
a[3] = -61271.89123310894 ;
a[4] = -24113.774471193552 ;
a[5] = 81251.17187865078 ;
a[6] = 2.6013991633662954 ;
a[7] = -47.47156475123484 ;
a[8] = 267.2414241330698 ;
a[9] = -1.0 ;
}
if (f==6){
a[0] = 55576.04523590207 ;
a[1] = 45585.559844881296 ;
a[2] = 736931.8212106228 ;
a[3] = 61411.1614453271 ;
a[4] = 24590.112847059965 ;
a[5] = 80733.09718613327 ;
a[6] = 4.488618038885761 ;
a[7] = 46.16711766895605 ;
a[8] = -266.8339533924591 ;
a[9] = -1.0 ;
}
if (f==7){
a[0] = 62606.78170695901 ;
a[1] = 51136.918440297246 ;
a[2] = 832667.3562214375 ;
a[3] = 63346.88473594934 ;
a[4] = -33744.77237933874 ;
a[5] = -72719.78481730819 ;
a[6] = -30.374276688322425 ;
a[7] = -24.862324570363853 ;
a[8] = -201.99640310788527 ;
a[9] = -1.0 ;
}
if (f==8){
a[0] = 55576.04523590207 ;
a[1] = 45585.559844881296 ;
a[2] = 736931.8212106228 ;
a[3] = -61411.1614453271 ;
a[4] = 24590.112847059965 ;
a[5] = -80733.09718613327 ;
a[6] = 4.488618038885761 ;
a[7] = -46.16711766895605 ;
a[8] = -266.8339533924591 ;
a[9] = -1.0 ;
}
}


//5th function = material

G4double G4BoundaryPhaseContrast::Mat(G4double a[],G4double i,G4double j,G4double k){
G4double fitted = a[0]* pow(i,2)+ a[1]*pow(j,2)+ a[2]*pow(k,2)+ a[3]*i*j +a[4]*i*k +a[5]*j*k + a[6]*i + a[7]*j + a[8]*k +a[9];
G4double mat;

if(fitted<=.0005){
mat= 0;
}
if(fitted>.0005){
mat =1;
}
return mat;
}




/////////////////////////////////////.........START CODE............./////////////////////////////////

G4BoundaryPhaseContrast::G4BoundaryPhaseContrast(const G4String &processName, G4ProcessType type) : G4VDiscreteProcess(processName, type) {

    if (verboseLevel > 0) {
        G4cout << GetProcessName() << " is created " << G4endl;
    }

    Material1 = NULL;
    Material2 = NULL;

   // kCarTolerance = G4GeometryTolerance::GetInstance()
     //      ->GetSurfaceTolerance();

    TotalMomentum = 0.;
    Rindex1 = Rindex2 = 1.;
    cost1 = cost2 = sint1 = sint2 = 0.;

}

G4BoundaryPhaseContrast::~G4BoundaryPhaseContrast() {}

G4VParticleChange *G4BoundaryPhaseContrast::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep) {


    aParticleChange.Initialize(aTrack);
    aParticleChange.ProposeVelocity(aTrack.GetVelocity());

    const G4Step *pStep = &aStep;

    const G4Step *hStep = G4ParallelWorldProcess::GetHyperStep();

    if (hStep) pStep = hStep;

///////////////////////////////////////////////////////////IS PARTICLE ON BOUNDARY?////////////////////////////


G4double xpoint =pStep->GetPostStepPoint()->GetPosition().x();
G4double ypoint =pStep->GetPostStepPoint()->GetPosition().y();
G4double zpoint =pStep->GetPostStepPoint()->GetPosition().z();

G4double f =G4BoundaryPhaseContrast::fit(xpoint,ypoint,zpoint);
G4double a[10]={0};
G4BoundaryPhaseContrast::index(f,a);



G4double marz = G4BoundaryPhaseContrast::boundary(a,xpoint,ypoint,zpoint);
    if (marz==1) {

//select material

            G4double ipoint =pStep->GetPreStepPoint()->GetPosition().x();
	    G4double jpoint =pStep->GetPreStepPoint()->GetPosition().y();
	    G4double kpoint =pStep->GetPreStepPoint()->GetPosition().z();

        G4double mat1=G4BoundaryPhaseContrast::Mat(a,ipoint,jpoint,kpoint);

        if(mat1==0){												//.............................................................//
           Material1= G4Material::GetMaterial("Liver");								//........materials should be defined in the macro file........//
														//.............................................................//
         }
        else if(mat1==1){
            Material1 =G4Material::GetMaterial("Vacuum");

        }
        G4double mat2=G4BoundaryPhaseContrast::Mat(a,xpoint,ypoint,zpoint);

        if(mat2==0){
            Material2= G4Material::GetMaterial("Liver");

        }
        else if(mat2==1){
            Material2 =G4Material::GetMaterial("Vacuum");

        }

    } else {

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    if (aTrack.GetStepLength() <= 1e-9 / 2) {
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();

    TotalMomentum     = aParticle->GetTotalMomentum();
    OldMomentum       = aParticle->GetMomentumDirection();
    OldPolarization   = aParticle->GetPolarization();




//find normal

    G4double difx=G4BoundaryPhaseContrast::diff(a,xpoint,ypoint,zpoint,1);
    G4double dify=G4BoundaryPhaseContrast::diff(a,xpoint,ypoint,zpoint,2);
    G4double difz=G4BoundaryPhaseContrast::diff(a,xpoint,ypoint,zpoint,3);
    G4double size=sqrt(pow(difx,2)+pow(dify,2)+pow(difz,2));
    G4double nx=difx/size;
    G4double ny=dify/size;
    G4double nz=difz/size;

   theGlobalNormal = G4ThreeVector(nx,ny,nz);

    Rindex1 = GetRindex(Material1, TotalMomentum);
    Rindex2 = GetRindex(Material2, TotalMomentum);

    G4double PdotN = OldMomentum * theGlobalNormal;

    cost1 = -PdotN;

    if (cost1 <1.0 - 1e-9) {        
        sint1 = std::sqrt(1. - cost1 * cost1);
        sint2 = sint1 * (Rindex1 / Rindex2);

    } else {
        sint1 = 0.0;
        sint2 = 0.0;
    }

    if (sint2 >= 1.0) {
        DoReflection();    
                       
    } else {
        if (cost1 > 0.0) {
            cost2 =  std::sqrt(1. - sint2 * sint2);
        } else {
            cost2 = -std::sqrt(1. - sint2 * sint2);
        }

        if (sint1 > .7 ) {               
            G4double alpha = cost1 - cost2 * (Rindex2 / Rindex1);
            NewMomentum = OldMomentum + alpha * theGlobalNormal;      
=
        } else {                                                            
            NewMomentum = OldMomentum;
            NewPolarization = OldPolarization;
        }
    }

    NewMomentum = NewMomentum.unit();
    NewPolarization = NewPolarization.unit();

    if ( verboseLevel > 0 ) {
        G4cout << " New Momentum Direction: " << NewMomentum     << G4endl;
        G4cout << " New Polarization:       " << NewPolarization << G4endl;
    }

    aParticleChange.ProposePolarization(OldPolarization);
    aParticleChange.ProposeMomentumDirection(NewMomentum);


    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

