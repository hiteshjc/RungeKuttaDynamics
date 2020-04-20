#include"global.h"
#include"search_for.h"
using namespace std;

void search_for(string str_in,string filename, string &str_ret, bool &found)
{
    // The input file has been formatted in a certain way 
    // We assume this form - For eg.
    // //////////////////////////////////////////////
    // hamiltonian=heisenberg
    // j_x=1.0
    // j_z=0.5
    ///...
    // //////////////////////////////////////////////
    // i.e. Format is
    // string_in=string_ret
    
    str_ret="";
    string str_read,str_bef,str_aft;
    
    found=false;
    int line=0;

    // Open file
    //cout<<"Defining ifstream"<<endl;
    ifstream fl(filename.c_str(),std::ios::in);
    //cout<<"Defined ifstream"<<endl;
    
    if (!fl)
    {   cerr<<"File "<<filename<<" could not be opened"<<endl;
        exit(1);
    }
    
    if (fl.is_open())
    {
        // Begin to read lines
        while ((! fl.eof()) and (not found))
        {
            getline(fl,str_read);
            line++;
            //cout<<"Line number I just read = "<<line<<endl;
            //cout<<"String I just read ="<<str_read<<endl;
            // Take input string and separator at "=" sign
            split(str_read,string("="),str_bef,str_aft);
            if (str_bef.compare(str_in)==0)
            {   
                fl.close();
                str_ret=str_aft;
                found=true;
            }
        }
        //if (not found) {cout<<" File "<<filename<<" was opened successfully but requested string ' "<< str_in<<" ' not found "<<endl;}
    }
}

void split(string &str, string sep_char, string &str_bef, string &str_aft)
{
    // Separate the string into two at the specified separating character
    size_t found;
    size_t length;
    int loc;

    char buffer[5000];
    // I am assuming the user wont put an entry greater than 100 charaters long! 
    // (50 was kinda not enough -given long name of files)

    found=str.find(sep_char);
    if (found!=string::npos)
    {
        loc=int(found);
        
        length=str.copy(buffer,loc,0);
        buffer[length]='\0';
        str_bef=buffer;
        
        length=str.copy(buffer,str.size()-loc,loc+1);
        buffer[length]='\0';
        str_aft=buffer;
    }
}
