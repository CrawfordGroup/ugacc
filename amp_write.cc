#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <psi4-dec.h>
#include "libparallel/ParallelPrinter.h"


/*
** Print the largest num_amps amplitudes to outfile.
*/
namespace psi { namespace ugacc {

#define AMP_CUTOFF 1e-8

struct onestack {
    double value;
    int i;
    int a;
};

struct twostack {
    double value;
    int i; int j;
    int a; int b;
};

void onestack_insert(struct onestack *stack, double value, int i, int a, 
    int level, int stacklen);
void twostack_insert(struct twostack *stack, double value, int i, int j, 
    int a, int b, int level, int stacklen);
void amp_write_T1(double **T1, int no, int nv, int length, std::string label, std::string OutFileRMR);
void amp_write_T2(double ****T2, int no, int nv, int length, std::string label, std::string OutFileRMR);

void amp_write(int num_amps, double **t1, double ****t2, int no, int nv, std::string name)
{
  std::string label;

  label = "\n\tLargest " + name + "1 Amplitudes:\n";
  amp_write_T1(t1, no, nv, num_amps, label, "outfile");
  label = "\n\tLargest " + name + "2 Amplitudes:\n";
  amp_write_T2(t2, no, nv, num_amps, label, "outfile");
}

void amp_write_T1(double **T1, int no, int nv, int length, std::string label, std::string out)
{
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            boost::shared_ptr<OutFile>(new OutFile(out)));

  int num2print=0;
  struct onestack *t1stack;

  t1stack = (struct onestack *) malloc(length * sizeof(struct onestack));
  for(int m=0; m < length; m++) { t1stack[m].value = 0; t1stack[m].i = 0; t1stack[m].a = 0; }

  int numt1 = no * nv;

  for(int i=0; i < no; i++) {
    for(int a=0; a < nv; a++) {
      double value = T1[i][a];
      for(int m=0; m < length; m++) {
        if((fabs(value) - fabs(t1stack[m].value)) > 1e-12) {
          onestack_insert(t1stack, value, i, a, m, length);
	  break;
	}
      }
    }
  }

  for(int m=0; m < ((numt1 < length) ? numt1 : length); m++)
    if(fabs(t1stack[m].value) > AMP_CUTOFF) num2print++;

  if(num2print) printer->Printf("%s", label.c_str());

  for(int m=0; m < ((numt1 < length) ? numt1 : length); m++)
    if(fabs(t1stack[m].value) > AMP_CUTOFF)
      printer->Printf("\t        %3d %3d %20.10f\n", t1stack[m].i, t1stack[m].a, t1stack[m].value);

  free(t1stack);
}

void onestack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen)
{
  int l;
  struct onestack temp;

  temp = stack[level];

  stack[level].value = value;
  stack[level].i = i;
  stack[level].a = a;

  value = temp.value;
  i = temp.i;
  a = temp.a;

  for(l=level; l < stacklen-1; l++) {
    temp = stack[l+1];

    stack[l+1].value = value;
    stack[l+1].i = i;
    stack[l+1].a = a;

    value = temp.value;
    i = temp.i;
    a = temp.a;
  }
}

void amp_write_T2(double ****T2, int no, int nv, int length, std::string label, std::string out)
{
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            boost::shared_ptr<OutFile>(new OutFile(out)));

  int num2print=0;
  struct twostack *t2stack;

  t2stack = (struct twostack *) malloc(length * sizeof(struct twostack));
  for(int m=0; m < length; m++) { 
    t2stack[m].value = 0; 
    t2stack[m].i = 0; t2stack[m].j = 0;
    t2stack[m].a = 0; t2stack[m].b = 0;
  }

  int numt2 = no * no * nv * nv;

  for(int i=0; i < no; i++) {
    for(int j=0; j < no; j++) {
      for(int a=0; a < nv; a++) {
        for(int b=0; b < nv; b++) {
          double value = T2[i][j][a][b];

          for(int m=0; m < length; m++) {
            if((fabs(value) - fabs(t2stack[m].value)) > 1e-19) {
	      twostack_insert(t2stack, value, i, j, a, b, m, length);
	      break;
	    }
	  }
        }
      }
    }
  }

  for(int m=0; m < ((numt2 < length) ? numt2 : length); m++)
    if(fabs(t2stack[m].value) > AMP_CUTOFF) num2print++;

  if(num2print) printer->Printf("%s", label.c_str());

  for(int m=0; m < ((numt2 < length) ? numt2 : length); m++)
    if(fabs(t2stack[m].value) > AMP_CUTOFF)
      printer->Printf("\t%3d %3d %3d %3d %20.10f\n", t2stack[m].i, t2stack[m].j, 
	      t2stack[m].a, t2stack[m].b, t2stack[m].value);

  free(t2stack);
}

void twostack_insert(struct twostack *stack, double value, int i, int j, int a, int b, 
		     int level, int stacklen)
{
  int l;
  struct twostack temp;

  temp = stack[level];

  stack[level].value = value;
  stack[level].i = i;
  stack[level].j = j;
  stack[level].a = a;
  stack[level].b = b;

  value = temp.value;
  i = temp.i;
  j = temp.j;
  a = temp.a;
  b = temp.b;

  for(l=level; l < stacklen-1; l++) {
    temp = stack[l+1];

    stack[l+1].value = value;
    stack[l+1].i = i;
    stack[l+1].j = j;
    stack[l+1].a = a;
    stack[l+1].b = b;

    value = temp.value;
    i = temp.i;
    j = temp.j;
    a = temp.a;
    b = temp.b;
  }
}

}} // namespace psi::ugacc
