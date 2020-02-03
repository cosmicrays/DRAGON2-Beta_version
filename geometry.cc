//
//  geometry.cc
//  DRAGON
//
//  Created by Luca Maccione on 17/07/13.
//
//

#include "geometry.h"
#include "grid.h"
#include "input.h"

TGeometry::TGeometry(TGrid* Coord_, Input* in) {
    Coord = Coord_;
    inp = in;
    vector<double> xvec = Coord->GetX();
    vector<double> yvec = Coord->GetY();
    vector<double> zvec = Coord->GetZ();
    
    density = vector<double>(xvec.size()*yvec.size()*zvec.size(), 1);
}

TSpiralGeometry::TSpiralGeometry(TGrid* Coord_, Input* in) : TGeometry(Coord_, in) {
  
    vector<double> xvec = Coord->GetX();
    vector<double> yvec = Coord->GetY();
    vector<double> zvec = Coord->GetZ();
  
    //density = vector<double>(xvec.size()*yvec.size()*zvec.size(), 0);
    
    for (int ix = 0; ix < xvec.size(); ix++) {
        double x1 = xvec[ix];
        for (int iy = 0; iy < yvec.size(); iy++) {
            double y1 = yvec[iy];
            for (int iz = 0; iz<zvec.size(); iz++) {
                double z1 = zvec[iz];
                
                double spiralarmdensity=0;
                
                //if(type == "None") spiralarmdensity=1.0; //MW130415; MW130513 take 1, not 0
                
                if(x1==0 && y1==0){ // very galactic center is not defined
                    x1=1.e-20;
                    y1=1.e-20;
                }
                
                if(in->SA_type == "CII" || in->SA_type == "NII")
                {
                    
                    // rotation counter clockwise: sun position now at x=8.3,y=0,z=0 (former position: x=0,y=8.3,z=0)
                    double x=x1*cos(-90.*Pi/180.)+y1*sin(-90.*Pi/180.);
                    double y=-1.*x1*sin(-90.*Pi/180.)+y1*cos(-90.*Pi/180.);
                    double z=z1;
                    
                    double x_sun=in->xobs*cos(-90.*Pi/180.)+in->yobs*sin(-90.*Pi/180.);
                    double y_sun=-1.*in->xobs*sin(-90.*Pi/180.)+in->yobs*cos(-90.*Pi/180.);
                    double z_sun=in->zobs;
                    
                    //constants
                    const int arms=4;
                    double sigma_r;
                    const double sigma_z=0.07;
                    double delta_i=15.*Pi/180.;
                    const double R3=2.9;
                    const double alpha[arms]={0.242,0.279,0.249,0.240};
                    const double a[arms]={0.246,0.608,0.449,0.378};
                    
                    double B[arms];
                    double B_max;
                    
                    if(in->SA_type=="CII"){
                        B[0] = 169.;
                        B[1] = 266.;
                        B[2] = 339.;
                        B[3] = 176.;
                        B_max=339.;
                    }
                    else if(in->SA_type=="NII"){
                        B[0] = 20.0;
                        B[1] = 31.4;
                        B[2] = 40.1;
                        B[3] = 20.8;
                        B_max=40.1;
                    }
                    else{ //MW130415
                        B[0] = 0.;
                        B[1] = 0.;
                        B[2] = 0.;
                        B[3] = 0.;
                        B_max= 0.;
                    }
                    
                    
                    double phi_i[arms];
                    double DELTA_phi[arms];
                    
                    double phi_i_sun[arms];
                    double DELTA_phi_sun[arms];
                    
                    double r;
                    double phi;
                    double x_tmp;
                    double y_tmp;
                    
                    double r_sun;
                    double phi_sun;
                    double x_sun_tmp;
                    double y_sun_tmp;
                    
                    
                    //calculate density
                    
                    r=sqrt(pow(x,2.0)+pow(y,2.0));
                    r_sun=sqrt(pow(x_sun,2.0)+pow(y_sun,2.0));
                    
                    phi=atan2(y,x);
                    phi_sun=atan2(y_sun,x_sun);
                    
                    if(phi<0) phi+=(2*Pi);
                    if(phi_sun<0) phi_sun+=(2*Pi);
                    
                    for(int i=0;i<arms;i++){
                        phi_i[i]=(log(r/a[i])/alpha[i]);
                        x_tmp=r*cos(phi_i[i]);
                        y_tmp=r*sin(phi_i[i]);
                        phi_i[i]=atan2(y_tmp,x_tmp);
                        if(phi_i[i]<0) phi_i[i]+=(2*Pi);
                        DELTA_phi[i]=fabs(phi-phi_i[i]);
                        if(DELTA_phi[i]>Pi) DELTA_phi[i]=2.*Pi-DELTA_phi[i];
                        
                        
                        phi_i_sun[i]=(log(r_sun/a[i])/alpha[i]);
                        x_sun_tmp=r_sun*cos(phi_i_sun[i]);
                        y_sun_tmp=r_sun*sin(phi_i_sun[i]);
                        phi_i_sun[i]=atan2(y_sun_tmp,x_sun_tmp);
                        if(phi_i_sun[i]<0) phi_i_sun[i]+=(2*Pi);
                        DELTA_phi_sun[i]=fabs(phi_sun-phi_i_sun[i]);
                        if(DELTA_phi_sun[i]>Pi) DELTA_phi_sun[i]=2.*Pi-DELTA_phi_sun[i];
                        
                    }
                    
                    if(in->SA_type=="CII"){
                        if(r<=R3) sigma_r=0.7;
                        else sigma_r=3.1;
                    }
                    if(in->SA_type=="NII"){
                        if(r<=R3) sigma_r=0.7;
                        else sigma_r=2.3;
                    }
                    
                    double res=0.;
                    double res_sun=0.;
                    
                    for(int i=0;i<arms;i++) res+=B[i]*exp(-1.*fabs(r-R3)/sigma_r)*exp(-1.*(pow(((DELTA_phi[i])/delta_i),2.0)))*exp(-1.*pow(z,2.0)/(2.*pow(sigma_z,2.0)));
                    for(int i=0;i<arms;i++) res_sun+=B[i]*exp(-1.*fabs(r_sun-R3)/sigma_r)*exp(-1.*(pow(((DELTA_phi_sun[i])/delta_i),2.0)))*exp(-1.*pow(z_sun,2.0)/(2.*pow(sigma_z,2.0)));
                    
                    res/=res_sun; //normalize to sun's position
                    
                    spiralarmdensity = res;
                }
                
                if(in->SA_type == "Daniele" || in->SA_type == "Blasi" || in->SA_type == "BlasiModel")
		//DG luglio2013 spiral arm implementation from Blasi papers
                {
                    //MW130623: concerning Simon's rotation (s.above) - do we need it here?
                    double x = x1;
                    double y = y1;
                    double zeta = z1;
                    
                    double radius = sqrt(x*x+y*y);
                    double theta    = atan2(y,x) ;
                    //
                    double spiral_arms = 0.;
	            int n_arms         = in->num_arms; 		 
                    double width       = in->spiral_width;
		    vector<double> arms_K  = in->arms_Kvec;          //vector<double>(n_arms,0.);
		    vector<double> arms_r0 = in->arms_r0vec;         //vector<double>(n_arms,0.);
		    vector<double> arms_theta0 = in->arms_theta0vec; //vector<dobule>(n_arms,0.);

		    double norm_gaussian = 1./(width*sqrt(2.*3.141));

		    for (int i_arms=0; i_arms < n_arms; i_arms++) {
			double arm_theta = (arms_K[i_arms] * log(radius/arms_r0[i_arms]) + arms_theta0[i_arms]);
		        int n1 = floor(fabs(arm_theta - theta)/(2.*M_PI));
                        int n2 =  ceil(fabs(arm_theta - theta)/(2.*M_PI));
                        double lower_r = arms_r0[i_arms] * exp((theta - arms_theta0[i_arms] + n1*2.*M_PI)/arms_K[i_arms]);
                        double upper_r = arms_r0[i_arms] * exp((theta - arms_theta0[i_arms] + n2*2.*M_PI)/arms_K[i_arms]);
                        double arm_distance = min(  fabs(radius-lower_r) , fabs(upper_r-radius) );
			spiral_arms += norm_gaussian * exp( -pow(arm_distance,2.)/(2.*width*width) );
	            }

                    
                    /*
                    double norma_K        = 4.25;
                    double norma_r0       = 3.48;
                    double norma_theta0   = 0.00;
                    double norma_theta    = (norma_K * log(radius/norma_r0) + norma_theta0);
                    double n1 = floor(fabs(norma_theta - theta)/(2.*M_PI));
                    double n2 = ceil(fabs(norma_theta - theta)/(2.*M_PI));
                    double lower_r = norma_r0 * exp((theta - norma_theta0 + n1*2.*M_PI)/norma_K);
                    double upper_r = norma_r0 * exp((theta - norma_theta0 + n2*2.*M_PI)/norma_K);
                    double norma_distance = min(  fabs(radius-lower_r) , fabs(upper_r-radius) );
                    //double norma_distance = sqrt( pow((x - radius*cos(norma_theta)), 2.) + pow((y - radius*sin(norma_theta)),2.) );
                    //
                    double carina_K       = 4.25;
                    double carina_r0      = 3.48;
                    double carina_theta0  = 3.14;
                    double carina_theta   = (carina_K * log(radius/carina_r0) + carina_theta0);
                    n1 = floor(fabs(carina_theta - theta)/(2.*M_PI));
                    n2 = ceil(fabs(carina_theta - theta)/(2.*M_PI));
                    lower_r = carina_r0 * exp((theta - carina_theta0 + n1*2.*M_PI)/carina_K);
                    upper_r = carina_r0 * exp((theta - carina_theta0 + n2*2.*M_PI)/carina_K);
                    double carina_distance = min(  fabs(radius-lower_r) , fabs(upper_r-radius) );
                    //double carina_distance = sqrt( pow((x - radius*cos(carina_theta)), 2.) + pow((y - radius*sin(carina_theta)),2.));
                    //
                    double perseus_K      = 4.89;
                    double perseus_r0     = 4.90;
                    double perseus_theta0 = 2.52;
                    double perseus_theta  = (perseus_K * log(radius/perseus_r0) + perseus_theta0);
                    n1 = floor(fabs(perseus_theta - theta)/(2.*M_PI));
                    n2 = ceil(fabs(perseus_theta - theta)/(2.*M_PI));
                    lower_r = perseus_r0 * exp((theta - perseus_theta0 + n1*2.*M_PI)/perseus_K);
                    upper_r = perseus_r0 * exp((theta - perseus_theta0 + n2*2.*M_PI)/perseus_K);
                    double perseus_distance = min(  fabs(radius-lower_r) , fabs(upper_r-radius) );
                    //double perseus_distance = sqrt( pow((x - radius*cos(perseus_theta)), 2.) + pow((y - radius*sin(perseus_theta)),2.) );
                    //
                    double scutum_K       = 4.89;
                    double scutum_r0      = 4.90;
                    double scutum_theta0  =-0.62;
                    double scutum_theta   = (scutum_K * log(radius/scutum_r0) + scutum_theta0);
                    n1 = floor(fabs(scutum_theta - theta)/(2.*M_PI));
                    n2 = ceil(fabs(scutum_theta - theta)/(2.*M_PI));
                    lower_r = scutum_r0 * exp((theta - scutum_theta0 + n1*2.*M_PI)/scutum_K);
                    upper_r = scutum_r0 * exp((theta - scutum_theta0 + n2*2.*M_PI)/scutum_K);
                    double scutum_distance = min(  fabs(radius-lower_r) , fabs(upper_r-radius) );
                    //double scutum_distance = sqrt( pow((x - radius*cos(scutum_theta)), 2.) + pow((y - radius*sin(scutum_theta)),2.) );
                    */
                    
                    spiralarmdensity = spiral_arms;
                }
                
                density[Coord->indexD(ix,iy,iz)] = spiralarmdensity;
            }
        }
    }
    return ;
}

double TGeometry::GetPattern(int x1, int y1,int z1){
    
    return density[Coord->indexD(x1,y1,z1)];
}


//UseSpiral:
void TGeometry::ApplySpiral(std::vector<double>& smoothdensity, double powerindex, double cut){
    
    double smooth_sum=0.;
    double Spiral_sum=0.;
    
    vector<double> Spiral_density(density.size(),0);
    for (int isp = 0; isp < Spiral_density.size(); isp++) {
        Spiral_density[isp] = pow(density[isp], powerindex) ;
        Spiral_density[isp] = min ( Spiral_density[isp], cut );
        Spiral_density[isp] = max ( Spiral_density[isp], 1./cut );
    }
    
    //weighting masses, not densities
    double volume;
    int ind;
    for(int ix=0; ix<Coord->GetDimX(); ix++)
        for(int iy=0; iy<Coord->GetDimY(); iy++)
            for(int iz=0; iz<Coord->GetDimZ(); iz++)
            {
                ind = Coord->indexD(ix, iy, iz);
                volume = Coord->GetDeltaX(ix) * Coord->GetDeltaY(iy) * Coord->GetDeltaZ(iz);
                
                smooth_sum += volume * smoothdensity[ind];
//MW130729: Choose between (A) and (B)
                Spiral_sum += volume * smoothdensity[ind] * Spiral_density[ind]; //(A)
//                 Spiral_sum += volume * Spiral_density[ind]; //(B)

// cout << "[MW-DEBUG SPIRAL] " << ix << " " << iy << " " << iz << " " << volume << " | " << smooth_sum << " " << Spiral_sum << endl;
            }
    
    // for(int i=0;i<density.size();i++){
    // gas_sum+=density[i];
    // Spiral_sum+=Spiral_density[i];
    // }
    
    if(Spiral_sum==0) return ;
    
    double boost=smooth_sum/Spiral_sum;
    
    for(int i=0;i<smoothdensity.size();i++){
        smoothdensity[i] *= boost*Spiral_density[i]; //(A)
//         smoothdensity[i] = boost*Spiral_density[i]; //(B)
        if(smoothdensity[i]<1.e-200){
            smoothdensity[i]=1.e-200;
        }
    }
    
    cout << "Spiral Arm Boost:  " << boost << endl;
    return ;
}
