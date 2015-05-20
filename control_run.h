/* esto es un comentario */

# define CHAIN 2 /* 0 - Without polymer chains
                    1 - With polymer chains (monolayers)
                    2 - Grafted polymer chains (brushes)
                    
                 */
# define CRITERIO 2 /* 1 -  formas de construir las monolayers (acerco el monomero más cercano a la pared)
                       2 -  Acerco la coordenadas x al monomero con la coorenada x más cercana a la pared derecha del poro)*/
# define pol 1 /* indicating polymer see units_adaptation.f90 */
#undef geometry  /* Select geometry 0 - Flat surface? 
# define geometry 1 /* Select geometry 0 - Flat surface?
                                     1 - Long Cylinder (1D)
                                     2 - Square section nanochannel 
                                     3 - Short 2D Cylinder 
                  */
# undef VDW 

