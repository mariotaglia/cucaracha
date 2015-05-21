/* esto es un comentario */

/* First Generalities about the system */
# define CHAIN 1 /* 0 - Without polymer chains
                    1 - With polymer chains (monolayers) 
                        This use pol == 1 See below
                    2 - Grafted polymer chains (brushes)
                 */

# define MUPOL 1 /* 0 - number of polymer regulated with sigma
                    1 - indicating use of mupol polymer see units_adaptation.f90 */

# define CRITERIO 2 /* 1 -  formas de construir las monolayers (acerco el monomero más cercano a la pared)
                       2 -  Acerco la coordenadas x al monomero con la coorenada x más cercana a la pared derecha del poro)*/


#undef geometry  /* Select geometry 0 - Flat surface? 
# define geometry 1 /* Select geometry 0 - Flat surface?
                                     1 - Long Cylinder (1D)
                                     2 - Square section nanochannel 
                                     3 - Short 2D Cylinder 
                  */
/* Second Details about the system and Interactions */
# define POL 0 /* Defines the type of polymer:
                      0 - PAH
                      1 - PMEP
                      2 - 
               */
# undef VDW 

