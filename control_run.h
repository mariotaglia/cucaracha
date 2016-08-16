/* esto es un comentario */

/* First Generalities about the system */
# define CHAIN 1 /* 0 - Without polymer chains
                    1 - Monolayer of polymer chains 
                        This use pol == 1 See below(?)
                    2 - Grafted polymer chains (brushes)
                 */

# define MUPOL 1 /* 0 - number of polymer is regulated with sigma (no polymer in bulk)
                    1 - indicating use of mupol polymer see units_adaptation.f90 */

# define CRITERIO 2 /* 1 -  Acerco el monomero más cercano a la pared. 
                            Busca Rmax y acerca a la pared respetando el theta (bueno para armar las figuras).
                       2 -  Acerco el monomero con la coordenada xmax 
                            a la pared derecha del poro. El y se deja en y=0. ( todas las conf. superpuestas)
                       3 -  Cadenas con Centro de Masa en las distintas layers (1 grado de libertad mas para las cadenas).
                            Esta opcion se usa con CHAIN==1
                       4 -  Cadenas con Centro de Masa en el centro del poro.
                            Esta opcion se usa con CHAIN==1
                    */
# define fsigmaq 1 /* 0 -  Sin regulacion de carga en la superficie
                      1 -  con equilibrio quimico (regulacion de carga) en la pared del poro
                   */

/* Second Details about the system and Interactions */
# define POL 0 /* Defines the type of polymer:
                      0 - PAH
                      1 - PMEP
                      2 - neutral 
               */
# undef PAHCL /* To use with POL == 0 (PAH), chemical eq. between PAH and Cl

/*Not implemented yet*/
# undef geometry  /* Select geometry 0 - Flat surface? 
# define geometry 1 /* Select geometry 0 - Flat surface?
                                     1 - Long Cylinder (1D)
                                     2 - Square section nanochannel 
                                     3 - Short 2D Cylinder 
                  */
/* END Not implemented */
# undef VDW 

/*Betamu constante*/
# undef BMu_const  /* Que hace esto? */

/**** For debugging change undef by define **** */
/* This lines prints information in standard output */

# undef fdebug /*Imprime cadenas*/
# define fdebug_set_pore_distrib
# undef fdebug_pxs

# undef fdebug_rota36
