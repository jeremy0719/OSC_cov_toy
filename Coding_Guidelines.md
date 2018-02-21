### Coding Guidelines      

This document describes the coding guidelines followed in the [PROSPECT oscillation fitting and sensitivity C++ package](https://github.com/PROSPECT-collaboration/OscSens_CovMatrix/tree/master/OscSensFitterCC) . It is suggested to follow these guidelines while implementing any future code.      

##### Why Guidelines?    

Having a set of guidelines will make it easier for future users to read and contribute to the code base. The main code base is in `C++11`, there are many tutorials and resources for users available on the web.      

##### Documentation    

Documentation is really important for the people who want to understand the code after you have written it. [Doxygen](http://www.stack.nl/~dimitri/doxygen/) has been used in order to create the documentation for the code.   
* At least a one-line description is included at the beginning of the class.     
* At least a one-line description is included for each data member.     
* Single line comments have been implemented by preceding comments with `///`.     
* Comment blocks have been placed between `/** ..... */`.         

Look at the [Doxygen Commenting Rules](http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html) to understand how to implement this in your code.

##### File Name Extensions     
Name extensions aren't required by current systems, but they will help people follow the structure of the code:     
* C++11 source files are named `filename.cc` with a single class defined in each file.  The file name matches the name of the class. If more than one class are defined in a file, the file name is matched with the most important class of the file.    
* C++11 header files are named `filename.hh`.      
 
##### Namespaces     
It can be convenient to include a `use` directive into your source code to save typing the namespace name.  However, these are not added to header files since they may cause unexpected namespace collisions.      

The standard library is explicitly referenced using the `std` namespace so that it is very clear where to find the documentation.  The `std` namespace are abbreviated using a `use` directive.      


##### Naming Conventions     
It can really help other people to understand your code if you use naming conventions.  For example, with naming conventions other people (including those who may be new to C++11) can tell at a glance if something is a function, a class or a variable.  Since the analysis is also heavily dependent on ROOT, these suggestions match the naming conventions established by ROOT.     

##### Type and Class Names    
A type name begins with a capital `T` (for example `TTree` or `TMath`).  In most cases, a single `T` was used, however in case of creation of a pure abstract base class or a pure mix-in class, `TV` or `TM` are to be used respectively.    

##### Class Field Names     
In order to distinguish local variables from class fields, it is easier to prefix every class field with an `f`.  If a field has a multi-word name, then each word is capitalized (e.g. fFieldWithLongName).    

##### Class Methods     
All public class methods start with capital letters. The package is written to be transparent and the variables and functions are made private only where there is no need for other classes to access the members and variables.   

##### Class Destructors    
Since you never know when a class might be used as a base class, always use a virtual destructor.    

##### Sub-Class Names    
A type name that is defined inside of a class should not begin with `T`, but should begin with a capital (e.g. `TMainObject::SubObject`). The exception to this rule is that typedefs which mimic STL types (e.g. iterators) should obey STL rules (e.g. `TMainObject::iterator`).    

##### Function, Variables and Parameter Names   

###### Function Names    
All member function names begins with a upper case.     

###### Local Variables      
Local variables and method parameters always start with a lower case letter. If a variable has a long name, then the following words are capitalized (e.g. aLongVariableName).    

###### Global Variables      
Global variables are avoided. If you really need that global variable use the prefix `g` and document it well.       

