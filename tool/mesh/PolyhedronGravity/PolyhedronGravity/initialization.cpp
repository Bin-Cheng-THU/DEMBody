//
//  initialization.cpp
//  
//
//  Created by Yang Yu on 10/25/17.
//
//

#include "initialization.h"

//

void LoadParas(char* parafile, string &polyfile, string &debrisfile, string &iniconfile, ODE_OPTION &odeopt, SOLAR &ss)
{
    int numchar, id, n;
    int flag;

    ifstream infile;
    string line, core, name, value;
    
    flag = 0;

    infile.open(parafile);
    
    while(getline(infile,line))
    {
        numchar = line.length();
        for (n=0; n<numchar; n++) if (line[n] == '#') break;
        core.assign(line,0,n);
        //
        core.erase(std::remove(core.begin(), core.end(), '\t'), core.end());
        core.erase(std::remove(core.begin(), core.end(), ' '), core.end());
        //
        if (! core.empty())
        {
            //
            numchar = core.length();
            id = core.find('=');
            assert((id>0) && (id+1<numchar));
            name.assign(core,0,id);
            n = numchar - id - 1;
            assert(n>0);
            value.assign(core,id+1,n);
            //
            if (name == "PolyhedronFile")       {polyfile = value;                              flag++;}
            else if (name == "DebrisFile")      {debrisfile = value;                            flag++;}
            else if (name == "IniConFile")      {iniconfile = value;                            flag++;}
            //
            else if (name == "FunOption")       {odeopt.FunOption = atof(value.c_str());        flag++;}
            else if (name == "StepSize")        {odeopt.StepSize = atof(value.c_str());         flag++;}
            else if (name == "EndStep")         {odeopt.EndStep = atof(value.c_str());         flag++;}
            else if (name == "StartStep")       {odeopt.StartStep = atof(value.c_str());        flag++;}
            else if (name == "OutputInterval")  {odeopt.OutputInterval = atof(value.c_str());   flag++;}
            //
            else if (name == "SolarTide")       {ss.SolarTide = atof(value.c_str());            flag++;}
            else if (name == "SolarPressure")   {ss.SolarPressure = atof(value.c_str());        flag++;}
            else if (name == "ReflectionRate")  {ss.Reflection = atof(value.c_str());           flag++;}
            //
            else {cout<< "Unrecoganized item "<<name<<"!"<<endl; exit(1);}
            //
        }
    }
    infile.close();
    
    if (flag != NumParameter)
    {
        cout<<flag<<endl;
        cout<<"Important parameter setting missing!"<<endl;
        exit(1);
    }
}

void LoadInicon(string iniconfile, DEBRIS &d, KEPORB &kep, PHASE &x, CLOUD &cl)
{
    int i, tmpi;
    const char* file;
    FILE * infile;
    
    CloudAlloc(d, cl);
    
    file = iniconfile.c_str();
    
    infile = fopen(file, "r");
    //
    fscanf(infile,"%lf%lf%lf%lf%lf%lf",&kep.SemiMajorAxis, &kep.Eccentricity, &kep.LongAscendNode, &kep.Inclination, &kep.ArgPeriapsis, &kep.MeanAnomaly);
    
    fscanf(infile,"%d%lf%lf%lf%lf%lf%lf%lf",&tmpi, &x.AsteroidOrien[0], &x.AsteroidOrien[1], &x.AsteroidOrien[2], \
                                                              &x.AsteroidOrien[3], &x.AsteroidAngVel[0], &x.AsteroidAngVel[1], &x.AsteroidAngVel[2]);
    
    for(i=0; i<d.NumDebris; i++) fscanf (infile,"%d%d%lf%lf%lf%lf%lf%lf",&tmpi, &cl.status[i], &cl.DebrisPos[i][0], &cl.DebrisPos[i][1], \
                                                                         &cl.DebrisPos[i][2], &cl.DebrisVel[i][0], &cl.DebrisVel[i][1], &cl.DebrisVel[i][2]);
    
    fclose (infile);
    
    QuatNorm(x.AsteroidOrien);
}


