/* esto es un comentario */

/* First Generalities about the system */
# define CHAIN 1 /* 0 - Without polymer chains
                    1 - Monolayer of polymer chains 
                        This use pol == 1 See below(?)
                    2 - Grafted polymer chains (brushes)
                 */

# define MUPOL 0 /* 0 - number of polymer is regulated with sigma (no polymer in bulk)
                    1 - indicating use of mupol polymer see units_adaptation.f90 */

# define CRITERIO 3 /* 1 -  Acerco el monomero más cercano a la pared)
                       2 -  Acerco la coordenadas x al monomero con la coorenada x 
                            más cercana a la pared derecha del poro) funciona para monolayers polimero adsorbido.
                       3 -  Cadenas con Centro de Masa en las distintas layers (1 grado de libertad mas para las cadenas).
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
/*Not implemented yet*/
#undef geometry  /* Select geometry 0 - Flat surface? 
# define geometry 1 /* Select geometry 0 - Flat surface?
                                     1 - Long Cylinder (1D)
                                     2 - Square section nanochannel 
                                     3 - Short 2D Cylinder 
                  */
/* END Not implemented */
# undef VDW 

/**** For debugging change undef by define **** */
/* This lines prints information in standard output */

# undef fdebug
# undef fdebug_set_pore_distrib
# undef fdebug_pxs

# undef fdebug_rota36
