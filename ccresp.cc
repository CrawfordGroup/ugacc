#include "ccpert.h"
#include "ccresp.h"

using namespace std;

namespace psi { namespace ugacc {


CCResp::CCResp(shared_ptr<CCPert> X, shared_ptr<CCPert> Y)
{

   H_ = X->H_;
   CC_ = X->CC_;
   HBAR_ = X->HBAR_;
   Lambda_ = X->CCLambda_;
   
   no_ = CC_->no_;
   nv_ = CC_->nv_;
   
   X_ = X;
   Y_ = Y;
}
CCResp::~CCResp()
{

}

double  CCResp::linresp(shared_ptr<CCPert> X, shared_ptr<CCPert> Y)
{

 int no = no_ ;
 int nv = nv_ ;

 double **l1 =   Lambda_->l1_;
 double ****l2 = Lambda_->l2_;

 double **X1_x = X->X1_ ;
 double **Y1_x = X->Y1_ ;

 double ****X2_x = X->X2_ ;
 double ****Y2_x = X->Y2_ ;

 double **X1_y = Y->X1_ ;
 double **Y1_y = Y->Y1_ ;

 double ****X2_y = Y->X2_ ;
 double ****Y2_y = Y->Y2_ ;

 double **pert_x = X->pert_;
 double **pert_y = Y->pert_;

 double **Aov_x = X->Aov_;
 double **Aoo_x = X->Aoo_;
 double **Avv_x = X->Avv_;
 double **Avo_x = X->Avo_;

 double ****Aovoo_x = X->Aovoo_;
 double ****Avvvo_x = X->Avvvo_;
 double ****Avvoo_x = X->Avvoo_;

 double **Aov_y = Y->Aov_;
 double **Aoo_y = Y->Aoo_;
 double **Avv_y = Y->Avv_;
 double **Avo_y = Y->Avo_;

 double ****Aovoo_y = Y->Aovoo_;
 double ****Avvvo_y = Y->Avvvo_;
 double ****Avvoo_y = Y->Avvoo_;

 double first=0;
 double second=0;

 double polar1=0, polar2=0, polar=0;


  for (int i=0; i< no; i++)
    for (int a=0; a< nv; a++){
      polar1 += Avo_x[a][i] * Y1_y[i][a] ;
      first += Avo_x[a][i] * Y1_y[i][a] ;
      for (int j=0; j< no; j++) 
        for (int b=0; b< nv; b++){
  	      polar1 += (0.50) * (Avvoo_x[a][b][i][j] + Avvoo_x[b][a][j][i])* Y2_y[i][j][a][b] ;
  	      second  += (0.50) * (Avvoo_x[a][b][i][j] + Avvoo_x[b][a][j][i])* Y2_y[i][j][a][b] ;
    } 
  }

  for (int i=0; i< no; i++)
    for (int a=0; a< nv; a++){
      for (int j=0; j< no; j++) 
        for (int b=0; b< nv; b++){
		  polar2 +=  l1[i][a] * pert_x[j][b+no] * (2.0 * X2_y[i][j][a][b] - X2_y[i][j][b][a]);  // 2. checked
        }
        for (int c=0; c< nv; c++)
		  polar2 +=  l1[i][a] * Avv_x[a][c] * X1_y[i][c];   // 2. checked
          for (int k=0; k< no; k++)
            polar2 -=  l1[i][a] * Aoo_x[k][i] * X1_y[k][a];   // 2. checked

        for (int j=0; j< no; j++)
          for (int b=0; b< nv; b++){         
            for (int c=0; c< nv; c++)
              polar2 += l2[i][j][b][c] * X1_y[i][a] * Avvvo_x[b][c][a][j] ;   // 3. checked
	
        for (int k=0; k< no; k++){
          polar2 -= 0.5 * l2[i][j][a][b] * X1_y[k][a] * Aovoo_x[k][b][i][j] ;    // 4. checked
          polar2 -= 0.5 * l2[i][j][a][b] * X1_y[k][b] * Aovoo_x[k][a][j][i] ;    // 4. checked
               
		  polar2 -= 0.5 * l2[i][j][a][b] * (X2_y[k][j][a][b]) * Aoo_x[k][i] ; // 4. checked
	      polar2 -= 0.5 * l2[i][j][a][b] * (X2_y[k][i][b][a]) * Aoo_x[k][j] ; // 4. checked 
        }

        for (int c=0; c< nv; c++){
          polar2 += 0.5 * l2[i][j][a][b] * (X2_y[i][j][a][c]) * Avv_x[b][c] ;  // 4. checked
          polar2 += 0.5 * l2[i][j][a][b] * (X2_y[j][i][b][c]) * Avv_x[a][c] ;  // 4. checked
        }
          } 
		
        polar2 += 2.0 * pert_x[i][a+no] * X1_y[i][a]; // 1. checked
   }


  outfile->Printf("\n polarizability first term: %20.14lf  polarizability second term: %20.14lf\n", polar1, polar2);
  //outfile->Printf("\n polarizability: %20.14lf\n", polar1 + polar2);
  return -1.0 * (polar1 + polar2) ;


 }

}} // psi::ugaccc
