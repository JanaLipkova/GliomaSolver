// reading a text file
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

int main () {
  
  vector<float> vec;
  float tmp;
  ifstream myfile ("input.txt");
  
  if (myfile.is_open()){
     //while (! myfile.eof() ) {
     while( myfile.good()){
     myfile >> tmp;
     if(myfile.eof()) break;
     vec.push_back(tmp);
     cout<<"tmp="<<tmp<<endl;
    }
  }

  for(vector<float>::iterator it = vec.begin(); it != vec.end(); ++it)
     printf("%f  \n", *it);
 
  cout<<"size="<<vec.size()<<std::endl;
  
  for(int i = 0; i<vec.size(); i++)
    printf("%f \n", vec[i]);
  
  return 0;
}
