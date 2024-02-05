
#include <Vc/Vc>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>

#include "easylogging++.h"
#include "semt/Semt.h"

#include "opendihu.h"

using namespace std;
/*
struct FiberPointData
{
  Vc::double_v values[8];
};
Vc_DECLARE_ALLOCATOR(FiberPointData)

int f(int a)
{
  double v = 0;
  for(long long int j = 0; j < 1e8; j++)
  {
    if (j % 17 == 5 || j % 13 == 2 || j%7 == 6)
      v += 2*a;
    if (j%2 == 0)
      v *= 2.;
    if (j%3 == 2)
      v /= 3.21;
  }
  return v;
}

class A
{
public:

  void init()
  {
    v_.resize(5);
    std::iota(v_.begin(), v_.end(), 0);
  }

  void print()
  {
    LOG(DEBUG) << "v: ";
    for (int i = 0; i < v_.size(); i++)
      LOG(DEBUG) << v_[i];
  }

  void getV(std::vector<int> &&v)
  {
    v = v_;
  }
protected:

  std::vector<int> v_;
};
*/
double gen() {
  double d = rand() * 1e-5;
  if (rand() % 2 == 0) {
    return -d;
  }
  return d;
}
int main(
    int argc,
    char *argv
        []) { /*
               LOG(DEBUG) << "debug" << f(1);
               LOG(INFO) << "info" << f(2);
               VLOG(1) << "vlog1 " << f(3);
               LOG(ERROR) << "logerror" << f(4);

               cout << "double vector size: " << Vc::double_v::size() << ",
               float vector size: " << Vc::float_v::size() << endl;

               // test Vc

               Vc::array<double,14> data;
               std::iota(data.begin(), data.end(), 10);

               LOG(INFO) << "data: " << data;
               //auto indexes = Vc::double_v::IndexType::IndexesFromZero();
               auto indexes = Vc::double_v([](int n) { return n; });
               Vc::double_v gathered = data[indexes];  // gathered == [0, 1, 2,
               ...] LOG(INFO) << "gathered: " << gathered;

               //Vc::double_v gathered2 = data[Vc::IndexesFromZero];  //
               gathered == [0, 1, 2, ...]
               //LOG(DEBUG) << "gathered2: " << gathered2;

               auto indexes2 = Vc::double_v([](int n) { return n+1; });
               LOG(INFO) << StringUtility::demangle(typeid(indexes2).name());
               Vc::double_v gathered3 = data[indexes2];  // gathered == [0, 1,
               2, ...] LOG(INFO) << "indexes2: " << indexes2 << ", gathered3: "
               << gathered3 << ",";

               for (int i = 0; i < int(ceil((float)data.size() /
               Vc::double_v::size())); i++)
               {
                 auto indexes = Vc::double_v([&i](int n) { return
               i*Vc::double_v::size() + n; }); Vc::double_v v = data[indexes];
                 LOG(INFO) << "indexes: " << indexes << ", v: " << v;
               }*/

  VecD<3, Vc::double_v> a{Vc::double_v(Vc::Zero), Vc::double_v(Vc::Zero),
                          Vc::double_v(Vc::One)};
  a[0][0] = 0.1;
  a[0][1] = 0.2;
  a[0][2] = 0.3;
  a[0][3] = 0.4;
  Vc::double_v inverseNorm = 1. / MathUtility::norm<3, Vc::double_v>(a);

  LOG(INFO) << a;
  LOG(INFO) << "inverseNorm: " << inverseNorm;
  LOG(INFO) << a * inverseNorm;
  VecD<3, Vc::double_v> normalizedA = a * inverseNorm;
  LOG(INFO) << MathUtility::norm<3, Vc::double_v>(normalizedA);

  MathUtility::normalize<3, Vc::double_v>(a);
  LOG(INFO) << a;
  exit(0);
#if 0  
  Vc::double_v v;
  for (int i = 0; i < Vc::double_v::size(); i++)
  {
    v[i] = i * (i%2==0? 1 : -1);
  }
  
  Vc::double_v epsilon = 0;
  LOG(INFO) << "epsilon: " << epsilon << ", v: " << v << ", abs: " << Vc::abs(v);
  Vc::where(Vc::abs(v) < 2.5) | epsilon = 1;
  LOG(INFO) << "epsilon: " << epsilon;

#if 0
  // measurement
  const long long N = 1e5;
  Vc::array<double,N> valuesA;
  Vc::array<double,N> valuesB;
  std::array<double,N> valuesArrayA;
  std::array<double,N> valuesArrayB;
  
  for (long long i = 0; i < N; i++)
  {
    valuesA[i] = (double)(gen());
    valuesB[i] = (double)(gen());
    valuesArrayA[i] = valuesA[i];
    valuesArrayB[i] = valuesB[i];
  }
  
  Vc::double_v cc = 0;
  
	auto start0 = chrono::steady_clock::now();

  const long long nEntries = int(ceil((float)valuesA.size() / Vc::double_v::size()));
  for (long long i = 0; i < nEntries; i++)
  {
    auto indexes = Vc::double_v([&i](int n) { return i*Vc::double_v::size() + n; });
    Vc::double_v a = valuesA[indexes];
    Vc::double_v b = valuesB[indexes];
    cc += a*b;
    //LOG(INFO) << "indexes: " << indexes << ", v: " << v;
  }
  
  double result = 0;
  for (long long i = 0; i < Vc::double_v::size(); i++)
  {
    result += cc[i];
  }
  
	auto end0 = chrono::steady_clock::now();
  
  LOG(INFO) << "resulta: " << result;
	auto start1 = chrono::steady_clock::now();
  
  result = 0;
  for (long long i = 0; i < valuesA.size(); i++)
  {
    result += valuesArrayA[i] * valuesArrayB[i];
  }
  LOG(INFO) << "resultb: " << result;
  
	auto end1 = chrono::steady_clock::now();

    
  LOG(INFO) << "duration Vc : " << chrono::duration_cast<chrono::microseconds>(end0 - start0).count();
  LOG(INFO) << "duration std: " << chrono::duration_cast<chrono::microseconds>(end1 - start1).count();
#endif
  
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpace;
  
  // create a function space
  std::string pythonConfig = "config = {}";
  DihuContext settings(argc, argv, pythonConfig);
  
  // create dummy node positions 
  std::vector<Vec3> nodePositions(5*7);
  std::array<element_no_t,2> nElementsPerCoordinateDirection{3,2}; 
  
  // create mesh by meshManager
  std::array<int,2> nRanks({1,1});  // this is completely serial
  
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace = DihuContext::meshManager()->template createFunctionSpace<DisplacementsFunctionSpace>("name", nodePositions, nElementsPerCoordinateDirection, nRanks);
  
  const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  const int nElementsLocal = displacementsFunctionSpace->nElementsLocal();
  
  
  // loop over elements
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += Vc::double_v::size())
  {
    Vc::int_v elementNoLocalv([elementNoLocal, nElementsLocal](int i)
    {
      return (i >= Vc::double_v::size() || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i); 
    });
  
    // get indices of element-local dofs
    
    std::array<Vc::int_v,nDisplacementsDofsPerElement> dofNosLocalv = displacementsFunctionSpace->getElementDofNosLocal(elementNoLocalv);
    
    LOG(INFO) << "element " << elementNoLocalv << ", type: " << StringUtility::demangle(typeid(decltype(elementNoLocalv)).name());
    LOG(INFO) << "dofNosLocal: " << dofNosLocalv;
  }

  /*
  Vc::Memory<Vc::double_v, 11> memory;
  std::iota(memory.begin(), memory.end(), 10);
  LOG(INFO) << memory;
  
  // scalar access:
  for (int i = 0; i < memory.entriesCount(); ++i) {
      double x = memory[i]; // read
      memory[i] = 2*x;     // write
  }
  // more explicit alternative:
  for (int i = 0; i < memory.entriesCount(); ++i) {
      double x = memory.scalar(i); // read
      memory.scalar(i) = 2*x;     // write
  }
  // vector access:
  for (int i = 0; i < memory.vectorsCount(); ++i) {
      Vc::double_v x = memory.vector(i); // read
      memory.vector(i) = 2*x;       // write
  }
  LOG(INFO) << memory;
  
  int nFiberPoints = 100;
  int nVcVectors = (nFiberPoints + Vc::double_v::Size - 1) / Vc::double_v::Size;

  std::vector<FiberPointData> values(nVcVectors);
  
  for (int i = 0; i < nVcVectors; i++)
  {
    values[i].values[0] = 8.0;
  }*/
  
/*  for (int i = 0; i < 256; i++)
  {
    std::cout << i << ":  " << char(i) << std::endl;
  }
*/
  //Vc::Memory<FiberPointData> data(nFiberPoints);
  
  // test SEMT 
  /*
  // Define an instance of variable with index 0 named x.
  DVAR(x, 0);

  //auto f = pow(x, INT(1)/INT(2));
  //auto f = ln(pow(x, INT(2)));
  auto f = pow(sqrt(x), INT(2));
  cout << "f = " << f << endl;
  
  auto df1 = deriv_t(f, x);
  cout << "f' = " << df1 << endl;

  auto g = x - (x + x);
  cout << "g = " << g << endl;
  
  auto dg1 = deriv_t(g, x);
  cout << "g' = " << dg1 << endl;

  */
  /*
  std::string pythonConfig = R"(
config = {
  "nElements": [2,2],
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  typedef FunctionSpace::FunctionSpace<
    Mesh::StructuredRegularFixedOfDimension<1>,
    BasisFunction::LagrangeOfOrder<1>
  > FunctionSpaceType;
  
  
  std::shared_ptr<FunctionSpaceType> functionSpace = settings.meshManager()->createFunctionSpace<FunctionSpaceType>("test", settings.getPythonConfig());
  
  const int nComponentsTest = 1;

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponentsTest>> test = functionSpace->template createFieldVariable<nComponentsTest>("test");
  
  //Vec globalVector = test->valuesGlobal();
  //Vec localVector = test->valuesLocal();
  
  //test->zeroEntries();
  //test->zeroGhostBuffer();
  test->setRepresentationGlobal();
  test->startGhostManipulation();
  LOG(DEBUG) << "inititalized, test: " << *test;


  int ownRankNo = functionSpace->meshPartition()->rankSubset()->ownRankNo();
  if (ownRankNo == 0)
  {
    std::array<int,2> indices({0,1});
    std::array<double,2> values({0.1,1.1});
      
    // here VecSetValues is used, not VecSetValuesLocal, because it is a plain vector
    test->setValues<2>(indices, values, INSERT_VALUES);
  }
  else
  {
    
    std::array<int,2> indices({0,1});
    std::array<double,2> values({0.2,1.2});
          
    test->setValues<2>(indices, values, INSERT_VALUES);
  }

  LOG(DEBUG) << "before finishGhostManipulation, test: " << *test;
  test->finishGhostManipulation();


  LOG(DEBUG) << "after finishGhostManipulation, test: " << *test;


  test->debug();*/
#endif
  return EXIT_SUCCESS;
}
