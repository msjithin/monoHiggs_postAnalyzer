#include <iostream>
#include <string>
#include "Burrito.cpp"

using namespace std;
class BuckysClass{
public:
  BuckysClass(string z){
    cout<<" Buckys class constructor..."<<endl;
    setName(z);
  }
  void setName(string x){
    name = x;
  }
  string getName(){
    return name;
  }
private: 
  string name;
  
};


int main()
{
  cout<<"Hello world!"<<endl;
  BuckysClass bo("Lucky ");
  cout<<bo.getName()<<endl;
  Burrito bro;
  bro.printText();
  return 0;
}


