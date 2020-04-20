# include "vector_utilities.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////

void extract_subset(std::vector<int> const &config, 
                    std::vector<int> const &touched_sites, 
                    std::vector<int> &sub_config)

{  
   sub_config.clear();
   std::vector<int>::const_iterator p;
   for (p=touched_sites.begin();p!=touched_sites.end();p++)
   {sub_config.push_back(config[(*p)]);}
}
///////////////////////////////////////////////////////////////////////////////////

bool check_subset(std::vector<int> const &vec1, 
                  std::vector<int> const &vec2)
{
   return includes(vec1.begin(),vec1.end(),vec2.begin(),vec2.end());  
}
///////////////////////////////////////////////////////////////////////////////////

static bool comparevecs(std::vector<int> const &vec1, 
                        std::vector<int> const &vec2)
{
     if (vec1.size()!=vec2.size()) {return (vec1.size()<vec2.size());}
     else
     {
        for (int i=0;i<vec1.size();i++)
        {
            if (vec1[i]!=vec2[i])
            {
                return (vec1[i]<vec2[i]);
            }
        }
     }
    return false;
}

///////////////////////////////////////////////////////////////////////////////////

void sort_vecs(std::vector< std::vector<int> > &lists)
{sort(lists.begin(),lists.end(),comparevecs);}

///////////////////////////////////////////////////////////////////////////////////

void remove_subsets(std::vector< std::vector<int> > &lists)
{
    std::vector< std::vector<int> >:: iterator p,q;
    bool exit;

    // Sort vectors
    sort_vecs(lists);
    // Then run a loop from first element to last checking subset condition
    p=lists.begin();
    while (p!=lists.end())
    {
        q=p+1;
        exit=false;
        while (q!=lists.end() and exit==false)
        {
            if (check_subset(*q,*p)) 
            {lists.erase(p);
             exit=true;
            }
            q++;
         }
        if (not exit)
        { p++;}
     }   
}

