#include<iostream>
#include<cctype>
#include<typeinfo>
#include<vector>
#include<climits>
#include<sstream>
#include "Convert_To_Double.hpp"
#include "Bad_Conversion.hpp"
#include "Common_Variables.hpp"

template <typename T> struct Arg_Obj{

  const std::string arg_name;
  std::string arg_name2;
  std::string arg_description;
  std::string arg_error_message;
  bool limit_error;
  T value;
  T upper_lim;
  T lower_lim;
  int count;
  Arg_Obj(std::string name, T v, T ll, T ul, std::string desc, std::string em, std::string& mp);
  Arg_Obj(std::string name, std::string name2, T v, T ll, T ul, std::string desc, std::string em, std::string& mp);
  void check_limits();
  void add_to_manpage(std::string& mp);
};

template <typename T> Arg_Obj<T>::Arg_Obj(std::string name, T v, T ll, T ul, std::string desc, std::string em, std::string& mp):arg_name(name),value(v),upper_lim(ul),lower_lim(ll),arg_description(desc),arg_error_message(em){
  arg_name2="";
  add_to_manpage(mp);
  count=0;
  check_limits();
}

template <typename T> Arg_Obj<T>::Arg_Obj(std::string name, std::string name2, T v, T ll, T ul, std::string desc, std::string em, std::string& mp):arg_name(name),arg_name2(name2),value(v),upper_lim(ul),lower_lim(ll),arg_description(desc),arg_error_message(em){
  add_to_manpage(mp);
  count=0;
  check_limits();
}

template <typename T> void Arg_Obj<T>::check_limits(){
  if((typeid(value) != typeid(std::string)) && (typeid(value) != typeid(char))){
    if(value < lower_lim || value > upper_lim){
      limit_error = true;
      std::cout<<"ERROR: "<<arg_error_message<<" Use --help or -h for more options.\n";
    } else{
      limit_error = false;
    }
  }
}

template <typename T> void Arg_Obj<T>::add_to_manpage(std::string& mp){

  std::string tidied_description;
  int N_nl_pos = arg_description.length()/65;
  std::stringstream ss_value;
  std::string str_value;
  ss_value << value;
  ss_value >> str_value;
  //int nl_pos = 0;
  
  tidied_description = "\t";
  tidied_description += arg_name;
  if(typeid(value) != typeid(bool))
    tidied_description += "x";
  if(arg_name2 != ""){
    tidied_description += ", ";
    tidied_description += arg_name2;
    tidied_description += "x";
  }
  tidied_description += "\n\t\t";

  if(arg_description.length() < 65){
    tidied_description += arg_description;
  }else{
    int start = 0;
    int delta = 64;
    int end = delta;
    for(int i=0;i<N_nl_pos;i++){

      if(arg_description[end] == ' '){
	tidied_description += arg_description.substr(start,delta);
	start = end+1;
	end += delta;
      }
      else{

	for(int j=delta; j > 0; j--){
	  if(arg_description[start+j] == ' '){
	    tidied_description += arg_description.substr(start,j);
	    start+=j+1;
	    end = start+delta;
	    break;
	  }
	}
      }
      tidied_description += "\n\t\t";
      if(end > arg_description.length())
	end = arg_description.length();
    }

    tidied_description += arg_description.substr(start,arg_description.length());
  }

  tidied_description += "\n\t\tDefault: ";
  tidied_description += str_value;
  tidied_description += "\n\n";

  mp += tidied_description;
}

