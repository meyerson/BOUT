#include <python2.6/Python.h>
#include <bout.hxx>
// #include <boutmain.hxx>

// #include <initialprofiles.hxx>
// #include <derivs.hxx>
// #include <interpolation.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int callPy(int argc, char *argv[])
{

  int rank, size,i;

  PyObject *pName, *pModule, *pDict, *pFunc,*pDir;
  PyObject *pArgs, *pValue;

  PyObject *sys_path; 
  PyObject *path,*path1, *path2, *path3; 

  
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
   
  // we can run serial code on the master node . . .
  if (rank == 0)
    {  

      printf("Running on processor %d \n",rank);
      Py_Initialize();
      
      PyRun_SimpleString("from time import time,ctime\n"
			 "print 'Today is',ctime(time())\n");

      FILE *fp = fopen("/home/cryosphere/BOUT/tools/pylib/py_try4.py","r+");
      
      sys_path = PySys_GetObject("path"); 
      if (sys_path == NULL) 
	return NULL; 
      path = PyString_FromString("/home/cryosphere/BOUT/tools/pylib/post_bout");
      if (path == NULL) 
	return NULL; 
      if (PyList_Append(sys_path, path) < 0) 
	return NULL; 
      Py_DECREF(path);
      
    
     
      PyRun_SimpleString("from time import time,ctime\n"
			 "print 'Today is',ctime(time())\n");
     

      pName = PyString_FromString(argv[0]); //module name
      
      output.write("pName: %s \n", argv[0]);
      output.write("path: %s \n", argv[2]);

      pModule = PyImport_ImportModule(argv[0]);
      PyObject *m_pDict = PyModule_GetDict(pModule); 
   
      Py_DECREF(pName);
      		  
      if (pModule != NULL) {
	
	pDir = PyObject_Dir(m_pDict);
	
	output.write(" PyList_Size(pDir): %i \n", PyList_Size(pDir));    
	
	pFunc = PyObject_GetAttrString(pModule,argv[1]); //single out a function from  a given module
	/* pFunc is a new reference */
	output.write("PyCallable_Check(pFunc): %i \n",PyCallable_Check(pFunc));
	
	pValue = PyString_FromString(argv[2]);
	pArgs = PyTuple_New(1);
	PyTuple_SetItem(pArgs, 0, pValue);
	
	if (pFunc && PyCallable_Check(pFunc)) {
	//   //parse the argument to pass to the python function
	   if (argc == 3){
	     if (argv[2] == NULL)
	       pArgs = NULL;
	   }
	   else if(argc ==2)
	     pArgs = NULL;
	   else {
	//     pArgs = PyTuple_New(argc - 2);
	//     for (i = 0; i < argc - 2; ++i) {
	//       //pValue = PyInt_FromLong(atoi(argv[i + 2]));
	//       pValue = PyString_FromString(argv[i + 2]);
	//       if (!pValue) {
	// 	Py_DECREF(pArgs);
	// 	Py_DECREF(pModule);
	//       fprintf(stderr, "Cannot convert argument\n");
	//       return 1;
	//       }
	//       PyTuple_SetItem(pArgs, i, pValue);
	//     }
	   }
	  
	  pValue = PyObject_CallObject(pFunc,pArgs);
  
	  if (pValue != NULL) {
	    printf("Result of call: %ld\n", PyInt_AsLong(pValue));
	    Py_XDECREF(pValue);
	  }
	  
	  else {
	    Py_XDECREF(pFunc);
	    Py_XDECREF(pModule);
	    PyErr_Print();
	    fprintf(stderr,"Call failed\n");
	    return 1;
	  }
	  Py_XDECREF(pFunc);
	  Py_XDECREF(pModule); 
	  
	}//close if function ok
	else {
	  if (PyErr_Occurred())
	    PyErr_Print();
	  fprintf(stderr, "Cannot find function \"%s\"\n",argv[1]);
        } //close if function not ok
      } //close if module ok
      else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n",argv[0]);
        return 1; 
      } //close if module not ok
      Py_Finalize();
    }// if cpu = 0
  

  
  return 0;
}
