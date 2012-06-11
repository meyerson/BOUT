/**************************************************************************
 *
 * Main BOUT++ functions
 * Adapted from the BOUT code by B.Dudson, University of York, Oct 2007
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

const char DEFAULT_DIR[] = "data";
const char DEFAULT_OPT[] = "BOUT.inp";
static char DEFAULT_GRID[] = "data/bout.grd.pdb";
static char help[] = "BOUT++: Uses finite difference methods to solve plasma fluid problems in curvilinear coordinates";


// MD5 Checksum passed at compile-time
#define CHECKSUM1_(x) #x
#define CHECKSUM_(x) CHECKSUM1_(x)
#define CHECKSUM CHECKSUM_(MD5SUM)

// Revision passed at compile time
#define REV1_(x) #x
#define REV_(x) REV1_(x)
#define REV REV_(REVISION)

#define GLOBALORIGIN

#include <bout.hxx>
#include <datafile.hxx>
#include <grid.hxx>
#include <solver.hxx>
#include <field2d.hxx>
#include <field3d.hxx>
#include <vector2d.hxx>
#include <vector3d.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <utils.hxx>
#include <invert_laplace.hxx>
#include <interpolation.hxx>
#include <boutmesh.hxx>
#include <boutexception.hxx>
#include <optionsreader.hxx>
#include <callPy.hxx>


#include <boundary_factory.hxx>
#include <boundary_standard.hxx>

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <string>
#include <list>
using std::string;
using std::list;

#ifdef BOUT_HAS_PETSC
#include <petsc.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef SIGHANDLE
#include <signal.h>
void bout_signal_handler(int sig);  // Handles segmentation faults
#endif

bool append = false;

BoutReal simtime;
int iteration;

const char *data_dir;

const string time_to_hms(BoutReal t);   // Converts to h:mm:ss.s format
char get_spin();                    // Produces a spinning bar

int bout_monitor(BoutReal t, int iter, int NOUT); // Function called by the solver each timestep


int bout_init(int argc, char **argv)
{
  int i, NOUT;
  BoutReal TIMESTEP;
  string grid_name;
  bool dump_float; // Output dump files as floats

  char dumpname[512];

  string grid_ext, dump_ext; ///< Extensions for restart and dump files

  //const char *data_dir; ///< Directory for data input/output
  const char *opt_file; ///< Filename for the options file
  
#ifdef CHECK
  int msg_point; ///< Used to return the message stack to a fixed point
#endif

#ifdef SIGHANDLE
  /// Set a signal handler for segmentation faults
  signal(SIGSEGV, bout_signal_handler);
#endif

  // Set default data directory
  data_dir = DEFAULT_DIR;
  opt_file = DEFAULT_OPT;

  /// Check command-line arguments
  /// NB: "restart" and "append" are now caught by options
  for (i=1;i<argc;i++) {
    if (strncasecmp(argv[i], "no", 2) == 0) Datafile::enabled = false; // No output. Used for scaling studies

    if (strncasecmp(argv[i], "-d", 2) == 0) {
      // Set data directory
      if (i+1 >= argc) {
        fprintf(stderr, "Useage is %s -d <data directory>\n", argv[0]);
        return 1;
      }
      i++;
      data_dir = argv[i];
    }
    if (strncasecmp(argv[i], "-f", 2) == 0) {
      // Set data directory
      if (i+1 >= argc) {
        fprintf(stderr, "Useage is %s -f <options filename>\n", argv[0]);
        return 1;
      }
      i++;
      opt_file = argv[i];
    }
  }

  // If no communicator is supplied, initialise MPI
  if (!BoutComm::getInstance()->isSet()) MPI_Init(&argc,&argv);

  // If BoutComm was set, then assume that MPI_Finalize is called elsewhere
  // but might need to revisit if that isn't the case

  int NPES, MYPE;
  MPI_Comm_size(BoutComm::get(), &NPES);
  MPI_Comm_rank(BoutComm::get(), &MYPE);

  /// Set up the output
  if (MYPE == 0) output.enable(); // Enable writing to stdout
  else output.disable(); // No writing to stdout

  /// Open an output file to echo everything to
  /// On processor 0 anything written to output will go to stdout and the file
  output.open("%s/BOUT.log.%d", data_dir, MYPE);

  /// Print intro
  output.write("\nBOUT++ version %.2f\n", BOUT_VERSION);
#ifdef REVISION
  output.write("Revision: %s\n", REV);
#endif
#ifdef MD5SUM
  output.write("MD5 checksum: %s\n", CHECKSUM);
#endif
  output.write("Code compiled on %s at %s\n\n", __DATE__, __TIME__);
  output.write("B.Dudson (University of York), M.Umansky (LLNL) 2007\n");
  output.write("Based on BOUT by Xueqiao Xu, 1999\n\n");

  output.write("Processor number: %d of %d\n\n", MYPE, NPES);

  /// Print compile-time options

  output.write("Compile-time options:\n");

#ifdef CHECK
  output.write("\tChecking enabled, level %d\n", CHECK);
#else
  output.write("\tChecking disabled\n");
#endif

#ifdef SIGHANDLE
  output.write("\tSignal handling enabled\n");
#else
  output.write("\tSignal handling disabled\n");
#endif

#ifdef PDBF
  output.write("\tPDB support enabled\n");
#else
  output.write("\tPDB support disabled\n");
#endif

#ifdef NCDF
  output.write("\tnetCDF support enabled\n");
#else
  output.write("\tnetCDF support disabled\n");
#endif

#ifdef _OPENMP
  output.write("\tOpenMP parallelisation enabled\n");
#else
  output.write("\tOpenMP parallelisation disabled\n");
#endif

#ifdef METRIC3D
  output.write("\tRUNNING IN 3D-METRIC MODE\n");
#endif

  /// Get the options tree
  Options *options = Options::getRoot();

  try {
      /// Load settings file
    OptionsReader *reader = OptionsReader::getInstance();
    reader->read(options, "%s/%s", data_dir, opt_file);

    // Get options override from command-line
    reader->parseCommandLine(options, argc, argv);
  }catch(BoutException *e) {
    output << "Error encountered during initialisation\n";
    output << e->what() << endl;
    return 1;
  }

  /// Start MPI
  //Dirty, dirty hack
#ifdef BOUT_HAS_PETSC
  SolverType type("");
  string solver_option;
  Options *solver_options = options->getSection("solver");
  solver_options->get("type", solver_option, "", false);
  if (!solver_option.empty()) type = solver_option.c_str();
  if (!(strcasecmp(type, SOLVERPETSC31) && strcasecmp(type, SOLVERPETSC))) PetscInitialize(&argc,&argv,PETSC_NULL,help);
#endif

  try {

    /////////////////////////////////////////////
    /// Get some settings

    OPTION(options, NOUT, 1);
    OPTION(options, TIMESTEP, 1.0);

    options->get("grid", grid_name, DEFAULT_GRID);

    OPTION(options, dump_float,   true);
    OPTION(options, non_uniform,  false);

    // Check if restarting
    bool restart;
    OPTION(options, restart, false);
    OPTION(options, append, false);

    /// Get file extensions
    options->get("dump_format", dump_ext, DEFAULT_FILE_EXT);

    /// Setup derivative methods
    if (derivs_init()) {
      output.write("Failed to initialise derivative methods. Aborting\n");
      return(1);
    }

    ////////////////////////////////////////////

    /// Create the mesh
    mesh = new BoutMesh();

    output.write("Setting grid format\n");

    /// Load the grid
    options->get("grid_format", grid_ext, "");
    if (grid_ext.empty()) mesh->addSource(new GridFile(data_format(grid_name.c_str()), grid_name.c_str())); // Guess format based on grid filename
    else mesh->addSource(new GridFile(data_format(grid_ext.c_str()), grid_name.c_str())); // User-specified format

    if (mesh->load()) {
      output << "Failed to read grid. Aborting\n";
      return 1;
    }

    /// Setup the boundaries
    BoundaryFactory* bndry = BoundaryFactory::getInstance();
    bndry->add(new BoundaryDirichlet(), "dirichlet");
    bndry->add(new BoundaryNeumann(), "neumann");
    bndry->add(new BoundaryRobin(), "robin");
    bndry->add(new BoundaryConstGradient(), "constgradient");
    bndry->add(new BoundaryZeroLaplace(), "zerolaplace");
    bndry->add(new BoundaryConstLaplace(), "constlaplace");
    bndry->addMod(new BoundaryRelax(10.), "relax");

    /// Set the file names
    sprintf(dumpname, "%s/BOUT.dmp.%d.%s", data_dir, MYPE, dump_ext.c_str());

    // Set file formats
    output.write("Setting file formats\n");
    dump.setFormat(data_format(dumpname));

    if (dump_float) dump.setLowPrecision(); // Down-convert to floats

    /// Add book-keeping variables to the output files

    // This is a temporary hack to get around datafile's limitations (fix soon)
    static BoutReal version = BOUT_VERSION;
    dump.add(version, "BOUT_VERSION", 0);
    dump.add(simtime, "t_array", 1); // Appends the time of dumps into an array
    dump.add(iteration, "iteration", 0);

    mesh->outputVars(dump);

    /// initialise Laplacian inversion code
    invert_init();

    output.write("Initialising physics module\n");
    /// Initialise physics module
#ifdef CHECK
    msg_point = msg_stack.push("Initialising physics module");
#endif


    /// Create the solver
    solver = Solver::Create();

    if (physics_init(restart)) {
      output.write("Failed to initialise physics. Aborting\n");
      return 1;
    }

#ifdef CHECK
    // Can't trust that the user won't leave messages on the stack
    msg_stack.pop(msg_point);
#endif

    /// Initialise the solver
    solver->setRestartDir(data_dir);
    if (solver->init(physics_run, argc, argv, restart, NOUT, TIMESTEP)) {
      output.write("Failed to initialise solver-> Aborting\n");
      return 1;
    }

    /// Set the filename for the dump files
    dump.setFilename(dumpname);

    if (!restart) {
      /// Write initial state as time-point 0

      // Run RHS once to ensure all variables set
      if (physics_run(0.0)) {
        output.write("Physics RHS call failed\n");
        return 1;
      }

      if (append) dump.append();
      else {
        dump.write();
        append = true;
      }
    }

  }catch(BoutException *e) {
    output << "Error encountered during initialisation\n";
    output << e->what() << endl;
    return 1;
  }

  return 0;
}

int bout_run()
{
  /// Run the solver
  output.write("Running simulation\n\n");
  int status,pbstatus;
  try {
    time_t start_time = time((time_t*) NULL);
    output.write("\nRun started at  : %s\n", ctime(&start_time));

    status = solver->run(bout_monitor);

    time_t end_time = time((time_t*) NULL);
    output.write("\nRun finished at  : %s\n", ctime(&end_time));
    output.write("Run time : ");
    
 
    
    

    int dt = end_time - start_time;
    int i = (int) (dt / (60.*60.));
    if (i > 0) {
      output.write("%d h ", i);
      dt -= i*60*60;
    }
    i = (int) (dt / 60.);
    if (i > 0) {
      output.write("%d m ", i);
      dt -= i*60;
    }
    output.write("%d s\n", dt);
  }catch(BoutException *e) {
    output << "Error encountered during initialisation\n";
    output << e->what() << endl;
    return 1;
  }

   
  try {
    /// Post-processing
    
    //output << "data_dir \n" << data_dir <<endl;
    //char *data_dir_cp;
    char *cpy_str = (char *)malloc(strlen(data_dir) + 1 * sizeof(char));
    //char *ptr;
    strcpy(cpy_str,data_dir);
    char* pbinput [3] = {"post_bout","save",cpy_str};
    ///char* pbinput [3] = {"post_bout","save",strcat(path_key,data_dir)};
    pbstatus = callPy(3,pbinput); 
  }catch(BoutException *e) {
    output << "Error encountered during post-processing\n";
    output << e->what() << endl;
    return 1;
  }
  
  return status;
}

int bout_finish()
{
  // Delete the solver
 

  delete solver;

  // Get and delete the mesh data sources
  list<GridDataSource*> source = mesh->getSources();
  for (list<GridDataSource*>::iterator it = source.begin(); it != source.end(); it++)
    delete *it;

  // Delete the mesh
  delete mesh;

  // Delete 3D field memory
  Field3D::cleanup();

  // Cleanup boundary factory
  BoundaryFactory::cleanup();

  // close MPI
  //Dirty, dirty hack
#ifdef BOUT_HAS_PETSC
  SolverType type("");
  string solver_option;
  Options *options = Options::getRoot();
  options = options->getSection("solver");
  options->get("type", solver_option, "", false);
 
 

  if (!solver_option.empty()) type = solver_option.c_str();

  if (!(strcasecmp(type, SOLVERPETSC31) && strcasecmp(type, SOLVERPETSC))) PetscFinalize();
  else if (!BoutComm::getInstance()->isSet()) MPI_Finalize();

  options = Options::getRoot();
#else
  // If BoutComm was set, then assume that MPI_Finalize is called elsewhere
  // but might need to revisit if that isn't the case
  if (!BoutComm::getInstance()->isSet()) MPI_Finalize();
#endif
  
 

  return 0;
}

/*!*************************************************************************
 * SOLUTION MONITOR FUNCTION
 *
 * Called each timestep by the solver
 **************************************************************************/

int bout_monitor(BoutReal t, int iter, int NOUT)
{
  // Data used for timing
  static bool first_time = true;
  static BoutReal wtime = 0.0;       ///< Wall-time since last output
  static BoutReal wall_limit, mpi_start_time; // Keep track of remaining wall time

#ifdef CHECK
  int msg_point = msg_stack.push("bout_monitor(%e, %d, %d)", t, iter, NOUT);
#endif

  // Set the global variables. This is done because they need to be
  // written to the output file before the first step (initial condition)
  simtime = t;
  iteration = iter;

  /// Write (append) dump file

  dump.write(append);
  append = true;

  /// Collect timing information
  int ncalls = solver->rhs_ncalls;
  BoutReal wtime_rhs   = solver->rhs_wtime;
  //BoutReal wtime_invert = 0.0; // wtime_invert is a global
  BoutReal wtime_comms = mesh->wtime_comms;  // Time spent communicating (part of RHS)
  BoutReal wtime_io    = Datafile::wtime;      // Time spend on I/O

  output.print("\r"); // Only goes to screen

  if (first_time) {
    /// First time the monitor has been called

    /// Get some options
    Options *options = Options::getRoot();
    OPTION(options, wall_limit, -1.0); // Wall time limit. By default, no limit
    wall_limit *= 60.0*60.0;  // Convert from hours to seconds

    /// Record the starting time
    mpi_start_time = MPI_Wtime(); // NB: Miss time for first step (can be big!)

    first_time = false;

    /// Print the column header for timing info
    output.write("Sim Time  |  RHS evals  | Wall Time |  Calc    Inv   Comm    I/O   SOLVER\n\n");

    /// Don't know wall time, so don't print timings
    output.write("%.3e      %5d        -     -    -    -    -    - \n", simtime, ncalls); // Everything else
  } else {
    wtime = MPI_Wtime() - wtime;

    output.write("%.3e      %5d       %.2e   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", 
        simtime, ncalls, wtime,
        100.0*(wtime_rhs - wtime_comms - wtime_invert)/wtime,
        100.*wtime_invert/wtime,  // Inversions
        100.0*wtime_comms/wtime,  // Communications
        100.*wtime_io/wtime,      // I/O
        100.*(wtime - wtime_io - wtime_rhs)/wtime); // Everything else
  }

  // This bit only to screen, not log file

  BoutReal t_elapsed = MPI_Wtime() - mpi_start_time;
  output.print("%c  Step %d of %d. Elapsed %s", get_spin(), iteration+1, NOUT, (time_to_hms(t_elapsed)).c_str());
  output.print(" ETA %s", (time_to_hms(wtime * ((BoutReal) (NOUT - iteration - 1)))).c_str());

  if (wall_limit > 0.0) {
    // Check if enough time left

    BoutReal t_remain = mpi_start_time + wall_limit - MPI_Wtime();
    if (t_remain < wtime) {
      // Less than 1 time-step left
      output.write("Only %e seconds left. Quitting\n", t_remain);

#ifdef CHECK
      msg_stack.pop(msg_point);
#endif
      return 1; // Return an error code to quit
    } else {
      output.print(" Wall %s", (time_to_hms(t_remain)).c_str());
    }
  }

  /// Reset clocks for next timestep

  mesh->wtime_comms = 0.0; // Reset communicator clock
  Datafile::wtime = 0.0;
  wtime_invert = 0.0;
  wtime = MPI_Wtime();

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}

/*!************************************************************************
 * Add variables to be solved
 **************************************************************************/

// NOTE: Here bout_solve is for backwards-compatibility. Eventually will be removed

void bout_solve(Field2D &var, const char *name)
{
  // Add to solver
  solver->add(var, ddt(var), name);
}

void bout_solve(Field3D &var, const char *name)
{
  solver->add(var, ddt(var), name);
}

void bout_solve(Vector2D &var, const char *name)
{
  solver->add(var, ddt(var), name);
  var.setBoundary(name);
}

void bout_solve(Vector3D &var, const char *name)
{
  solver->add(var, ddt(var), name);
}

/*!************************************************************************
 * Add constraints
 **************************************************************************/

bool bout_constrain(Field3D &var, Field3D &F_var, const char *name)
{
  if (!solver->constraints()) return false; // Doesn't support constraints

  // Add to solver
  solver->constraint(var, F_var, name);

  return true;
}

/**************************************************************************
 * Global error handling
 **************************************************************************/

/// Print an error message and exit
void bout_error()
{
  bout_error(NULL);
}

void bout_error(const char *str)
{
  output.write("****** ERROR CAUGHT ******\n");

  if (str != NULL) output.write(str);

  output.write("\n");

#ifdef CHECK
  msg_stack.dump();
#endif

  MPI_Abort(BoutComm::get(), 1);

  exit(1);
}

#ifdef SIGHANDLE
/// Signal handler - catch segfaults
void bout_signal_handler(int sig)
{
  /// Set signal handler back to default to prevent possible infinite loop
  signal(SIGSEGV, SIG_DFL);

  output.write("\n****** SEGMENTATION FAULT CAUGHT ******\n\n");
#ifdef CHECK
  /// Print out the message stack to help debugging
  msg_stack.dump();
#else
  output.write("Enable checking (-DCHECK flag) to get a trace\n");
#endif

  exit(sig);
}
#endif

/**************************************************************************
 * Utilities
 **************************************************************************/

/// Write a time in h:mm:ss.s format
const string time_to_hms(BoutReal t)
{
  int h, m;

  h = (int) (t / 3600); t -= 3600.*((BoutReal) h);
  m = (int) (t / 60);   t -= 60 * ((BoutReal) m);

  char buffer[256];
  sprintf(buffer,"%d:%02d:%04.1f", h, m, t);

  return string(buffer);
}

/// Produce a spinning bar character
char get_spin()
{
  static int i = 0;
  char c = '|'; // Doesn't need to be assigned; squash warning

  switch(i) {
  case 0: 
    c = '|'; break;
  case 1:
    c = '/'; break;
  case 2:
    c = '-'; break;
  case 3:
    c = '\\'; break;
  }
  i = (i+1) % 4;
  return c;
}
