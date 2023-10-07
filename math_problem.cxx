// To run: root math_problem.cxx
// 
// Implementation of simple physics problem:
//   mass 1, originally at rest, v1_0=0
//   mass 2, originally at speed -v2_0 (from right to left)
// Placed as:
// 
// ||| 
// |||____|m1|____<--|m2|___
// 
// ||| = a wall that reflects m1 without energy loss
//
// Statement of the problem: the number of collisions needed to reflect back m2
// tends asymtotically to pi/2*sqrt(m2)
//
// I have implemented the following:
//   --Global variables: masses and initial speed of m1 and m2
//                       mech_energy = total mechanical energy, conserved
//                           it must be updated if m1,m2, v1_0, v2_0 are changed!
//			 tolerance = tolerance for the convergence
//   --calculate_step(...): calculate the change of speed of m1 and m2, given initial speeds
//   --count_reflections(): inside a while loop invoke calculate_step and reflects v1 
//			   (as it would be reflected back from the wall). 
//			   The loop breaks if the final speed of m2 changes less than "tolerance"
//   --math_problem(): runs the count_reflections() for different masses m2
//
//   Analytical solution: found later in youtube, https://youtu.be/jsYwFizhncE
//			  tldr, solve the problem in the phase space


double m1 = 1;
double m2 = 1;
double v1_0 = 0;
double v2_0 = -100000; /*from right to left, negative sign!*/
double mech_energy = 0.5*m1*v1_0*v1_0 + 0.5*m2*v2_0*v2_0; 

void update_energy(){mech_energy = 0.5*m1*v1_0*v1_0 + 0.5*m2*v2_0*v2_0; }

double tolerance = 1e-15;

bool MYDEBUG = false;

void calculate_step( double v1_pre, double v2_pre, double & v1_post, double & v2_post )
{

if(MYDEBUG) std::cout << "   *   *   *   *   *   *" << std::endl; 
if(MYDEBUG) std::cout << "v1_pre: " << v1_pre << std::endl; 
if(MYDEBUG) std::cout << "v2_pre: " << v2_pre << std::endl; 

double w = m1*v1_pre+m2*v2_pre;
double a = 0.5*m1*(1+m1/m2);
double b= -w*m1/m2;
double c = 0.5/m2*w*w - mech_energy;

double discr = b*b-4*a*c;

if(MYDEBUG) std::cout << "a: " << a << std::endl; 
if(MYDEBUG) std::cout << "b: " << b << std::endl; 
if(MYDEBUG) std::cout << "c: " << c << std::endl; 
if(MYDEBUG) std::cout << "discr: " << discr << std::endl; 

// x = v1_post, 2 solutions
double x_plus  = (-b + sqrt(b*b-4*a*c))/(2*a);
double x_minus = (-b - sqrt(b*b-4*a*c))/(2*a);

if(MYDEBUG) std::cout << "x_plus: " << x_plus << std::endl; 
if(MYDEBUG) std::cout << "x_minus: " << x_minus << std::endl; 

//after some test, it seems the good solution is x_minus
v1_post = x_minus;
v2_post = (w -m1*v1_post) / m2;

if(MYDEBUG) std::cout << "v1_post: " << v1_post << std::endl; 
if(MYDEBUG) std::cout << "v2_post: " << v2_post << std::endl; 
if(MYDEBUG) std::cin.ignore();
}


int count_reflections()
{
//set initial conditions
double v1_pre = 0;
double v2_pre = v2_0;
double v1_post = 0;
double v2_post = 0;

int n_reflections = 0;

while(true)
{
++n_reflections;
calculate_step( v1_pre, v2_pre, v1_post, v2_post );
v1_pre = -v1_post;
if( v2_post > 0 && fabs(v2_pre - v2_post) < tolerance ) break;
v2_pre = v2_post;
v1_post = 0;
v2_post = 0;
}

return n_reflections;
}

void math_problem()
{


std::vector<double> m_vector = {1,1e2,1e4,1e6,1e8,1e12,1e14};
std::vector<double> n_vector;
std::vector<double> time_vector;

for( auto m : m_vector)
{
std::cout << "Starting " << m << "...\n";
m2 = m;
// the final result actually do not depend on v2_0
// but convergence is faster if v2_0 is larger for larger m2 values
v2_0 = -10*sqrt(m); 
update_energy();

// monitor CPU/real time
TStopwatch a;
a.Start();
int ncollisions = count_reflections();
a.Stop();

time_vector.push_back(a.RealTime());
n_vector.push_back( ncollisions );
// std::cout << count_reflections() << std::endl;

std::cout << "Starting " << m << "... Done!\n";
}

std::cout << "    *    *    *    *    *    *    *    *\n" << std::endl; 
std::cout << "  m2" << "\t\t" << "  2*collisions/sqrt(m2)\t\t Real Time (s)" << std::endl;


for( int i = 0; i <m_vector.size(); ++i)
{
	std::cout <<  std::scientific << std::setprecision(2) << m_vector[i] << "\t";
	std::cout << std::fixed << std::setprecision(7) <<2*n_vector[i]/sqrt(m_vector[i])<< "\t\t\t";
	std::cout << std::fixed << time_vector[i] << std::endl;
}
}
