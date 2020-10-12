#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

int main ()
{
  std::ifstream file1("theirs_rootfile/eventAnalysis.txt");
  std::ifstream file2("mine_rootfile/dyjets.txt");
  std::string line1;
  std::string line2;
  std::vector<std::string> myevents1;
  std::vector<std::string> myevents2;
  std::vector<std::string> events1;
  std::vector<std::string> events2;
  while (std::getline(file2, line2)) {
    //std::cout << line << "\n";
    events2.clear();
    std::istringstream iss2(line2);
    std::string token2;
    while(std::getline(iss2, token2, '\t'))   // but we can specify a different one
      events2.push_back(token2);
    //std::cout<<"token = "<<events2[0]<<std::endl;
    myevents2.push_back(events2[0]);
  }
  
  while (std::getline(file1, line1)) {
    //std::cout << line << "\n";
    events1.clear();
    std::istringstream iss1(line1);
    std::string token1;
    while(std::getline(iss1, token1, '\t'))   // but we can specify a different one
      events1.push_back(token1);
    myevents1.push_back(events1[0]);
  }
  
  for(int i=0; i<myevents2.size(); i++){
    if(std::find(myevents1.begin(), myevents1.end(), myevents2[i]) != myevents1.end()) {
      /* v contains x */
    } else {
      std::cout<<"token = "<<myevents2[i]<<std::endl;
    }
  }

}
  


