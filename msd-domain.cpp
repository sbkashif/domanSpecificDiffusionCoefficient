/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iterator>
#include <algorithm>


#include <gromacs/trajectoryanalysis.h>

using namespace gmx;

/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
class AnalysisTemplate : public TrajectoryAnalysisModule
{
    public:
        AnalysisTemplate();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();
        // Adding new variables
        std::vector<double>                 frame_number;
        std::vector<double>                 frame_time;
        std::vector<double>                 msdt;
        std::vector<std::vector<double>>                 coord_sol_x;
        std::vector<std::vector<double>>                 coord_sol_y;
        std::vector<std::vector<double>>                 coord_sol_z;
        
        double                          zmax_i=10, zmin_i=-10;
        double                          zmax,zmin,comz;
        double                          boxx,boxy,boxz;
        double                          dis[4],pos_x,pos_y,pos_z;
        double**                        msd_pos;
        double**                        init_pos;
        double*                         msd_time;
        int*                            msd_index;
        int                             atom_num;
        double                          d_time = 0;
        int                             d_frame = 0;
        int                             timeflag=1;
        double                          fr_time_start;
        int                             n_start=10;

        std::vector<double>                     dym_ele{std::vector<double>(8,0)};
        std::vector<std::vector<std::vector<std::vector<double> > > >       record_dym;
        std::vector<std::vector<double> >                                   record_dym_i;
        std::vector<std::vector<std::vector<double> > >                     runoutofnames;
        std::vector<std::vector<double> >                                   dym_i;                              //write outp    ut
        std::vector<double>                                                 dym_i1{std::vector<double>(8,0)} ;  //write outp    ut
        std::vector<int>                                                    count;
        int nf=0;
        int ille_a=0;



    private:
        class ModuleData;

        std::string                      fnMSD_="msd_code";
        double                           cutoff_;
        Selection                        basesel_;
        Selection                        sel_;

        AnalysisNeighborhood             nb_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;
        // Adding the pointer about data points
        AnalysisDataPlotSettings         plotSettings_;
};


AnalysisTemplate::AnalysisTemplate()
    : cutoff_(0.0)
{
    registerAnalysisDataset(&data_, "avedist");
}


void
AnalysisTemplate::initOptions(IOptionsContainer          *options,
                              TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "This is a template for writing your own analysis tools for",
        "GROMACS. The advantage of using GROMACS for this is that you",
        "have access to all information in the topology, and your",
        "program will be able to handle all types of coordinates and",
        "trajectory files supported by GROMACS. In addition,",
        "you get a lot of functionality for free from the trajectory",
        "analysis library, including support for flexible dynamic",
        "selections. Go ahead an try it![PAR]",
        "To get started with implementing your own analysis program,",
        "follow the instructions in the README file provided.",
        "This template implements a simple analysis programs that calculates",
        "average distances from a reference group to one or more",
        "analysis groups."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("o")
                           .filetype(eftPlot).outputFile()
                           .store(&fnMSD_).defaultBasename("avedist")
                           .description("Average distances from reference group"));

    options->addOption(SelectionOption("Base")
                           .store(&basesel_).required()
                           .description("Base group to calculate distances from"));

    options->addOption(SelectionOption("select")
                           .store(&sel_).required()
                           .description("Groups to calculate MSD(OW) for"));

    options->addOption(DoubleOption("cutoff").store(&cutoff_)
                           .description("Cutoff for distance calculation (0 = no cutoff)"));

    options->addOption(DoubleOption("maxz").store(&zmax_i)
                           .description("Max z for plane with respect to the base"));

    options->addOption(DoubleOption("minz").store(&zmin_i)
                           .description("Min z for plane with respect to the base"));

    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);

}


void
AnalysisTemplate::initAnalysis(const TrajectoryAnalysisSettings &settings,
                               const TopologyInformation         & /*top*/)
{
    nb_.setCutoff(cutoff_);

    data_.setColumnCount(0, 1);

    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);

    std::cout << "You have picked:"<<basesel_.name()<<"as your base"<<std::endl;
    std::cout<< "You have picked:"<<sel_.name()<<"as your group for MSD calculation"<<std::endl;    
    
    if (!fnMSD_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnMSD_);
        plotm->setTitle("MSD vs Time");
        plotm->setXAxisIsTime();
        plotm->setYLabel("MSD");
        data_.addModule(plotm);
    }

    atom_num=sel_.atomCount();
    msd_index=new int [atom_num]();

    for (int i=0; i<atom_num;i++)
    {
        msd_index[i]=0;
        record_dym.push_back(runoutofnames);
    }
    std::cout << "Total number of atoms processed is:" << atom_num <<std::endl;
    count.assign(atom_num,-1);
    std::cout<<count[atom_num-1]<<std::endl;
}


void
AnalysisTemplate::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle         dh     = pdata->dataHandle(data_);
    const Selection           &sel = pdata->parallelSelection(sel_);
    const Selection           &basesel = pdata->parallelSelection(basesel_);

    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, sel);
    //dh.startFrame(frnr, fr.time);
    int nr_base=basesel.posCount();
    int nr_sel=sel.posCount();
    std::vector<double> coord_frame_x;
    std::vector<double> coord_frame_y;
    std::vector<double> coord_frame_z;

    if(timeflag==1)
        {
        timeflag=0;
        fr_time_start=fr.time;
        }       

    dym_i1[0]=fr.time-fr_time_start;
    dym_i.push_back(dym_i1);
    boxx=fr.box[0][0];
    boxy=fr.box[1][1];
    boxz=fr.box[2][2];
    comz = 0;



        for (int i = 0; i < nr_base; ++i)
        {
            SelectionPosition p = basesel.position(i);
            pos_z=p.x()[2];
            
            //std::cout<<"Atom:"<<i<<" Position:"<<pos_z<<std::endl;
            if(pos_z<0 || pos_z>=boxz)
            {
                pos_z=pos_z-boxz*floor(pos_z/boxz);
            }
            comz+=pos_z;
        }

        //std::cout<<" Before:"<<comz<<std::endl;
        comz /=nr_base; 
        //std::cout<<" After:"<<comz<<std::endl;
        zmax=comz+zmax_i;
        zmin=comz+zmin_i;
        //std::cout<<"Zmin:"<<zmin<<" Zmax:"<<zmax<<std::endl;
        
    for(int i=0;i<nr_sel;i++)
    {
        SelectionPosition p = sel.position(i);
        pos_x=p.x()[0];
        pos_y=p.x()[1];
        pos_z=p.x()[2];
        
        dym_ele[0]=frnr;
        dym_ele[1]=fr.time;
        
        dym_ele[2]=p.x()[0];
        dym_ele[3]=p.x()[1];
        dym_ele[4]=p.x()[2];

        if(pos_z<0 || pos_z>=boxz)
        {
            pos_z=pos_z-boxz*floor(pos_z/boxz);
        }

        if (pos_z<=zmax && pos_z>=zmin)
        {   
            coord_frame_x.push_back(dym_ele[2]);
            coord_frame_y.push_back(dym_ele[3]);
            coord_frame_z.push_back(dym_ele[4]);
        }
        else
        {
            ille_a++;
            coord_frame_x.push_back(-1000.0);
            coord_frame_y.push_back(-1000.0);
            coord_frame_z.push_back(-1000.0);
        }     

    }
    
    coord_sol_x.push_back(coord_frame_x);
    coord_sol_y.push_back(coord_frame_y);
    coord_sol_z.push_back(coord_frame_z);
    frame_time.push_back(fr.time);
        //frave /= nr;
        //dh.setPoint(g, frave);
    }
    //dh.finishFrame();
void
AnalysisTemplate::finishAnalysis(int nframes)
{
    std::vector<double>                 msd_x(nframes,0.0);
    std::vector<double>                 msd_y(nframes,0.0);
    std::vector<double>                 msd_z(nframes,0.0);
    std::vector<double>                 norm(nframes,0.0);
    
    std::cout << "Size:" << coord_sol_x[0].size() << std::endl;
    std::cout << "Size2:" << coord_sol_y[1].size() << std::endl;
    
    for(int i=0;i<nframes;i++)
    {
        //if (i%int(nframes/10)==0)
        //{
           std::cerr<<i<<" out of"<<nframes<<" frames processed"<< std::endl;
        //}
        for (int j=1;j<(nframes-i);j++)
        {
            for (int k=0;k<atom_num;k++)    
            {
                if (coord_sol_z[j+i][k] != -1000.0 && coord_sol_z[i][k]!= -1000.0)
                {
                    msd_x[j] += (coord_sol_x[j+i][k]-coord_sol_x[i][k])*(coord_sol_x[j+i][k]-coord_sol_x[i][k]); 
                    msd_y[j] += (coord_sol_y[j+i][k]-coord_sol_y[i][k])*(coord_sol_y[j+i][k]-coord_sol_y[i][k]); 
                    msd_z[j] += (coord_sol_z[j+i][k]-coord_sol_z[i][k])*(coord_sol_z[j+i][k]-coord_sol_z[i][k]);
                norm[j]++;
                }
            }
        }
    }

    norm[0]=1; 
    for(int k = 0; k < nframes; ++k)
        {
            std::cout << "Norm:" << norm[k] <<std::endl;
            msd_x[k] /= norm[k];
            msd_y[k] /= norm[k];
            msd_z[k] /= norm[k];
            msdt.push_back(msd_x[k] + msd_y[k] + msd_z[k]);
        }
}


void
AnalysisTemplate::writeOutput()
{
    // We print out the average of the mean distances for each group.
    std::string filename2;
    filename2 =fnMSD_.c_str();
    FILE * myfile2;
    myfile2 = fopen(filename2.c_str(), "a");
    fprintf (myfile2,"%s\n", "@ time difference     MSD_x       MSD_y     MSD_xyz" );
    for (int i=0; i<msdt.size(); i++)
    {
        fprintf (myfile2, "%12.3f%14.5f\n",frame_time[i],msdt[i]);
    }
    fclose(myfile2);


}

/*! \brief
 * The main function for the analysis template.
 */
int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AnalysisTemplate>(argc, argv);
}
