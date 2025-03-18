def calculate_hbonds(universe, donors, acceptors, dist=3.5, angle=150, start=None, stop=None, step=None):
    """
    Analiza wiązań wodorowych: liczymy ile jest H-bonds w czasie.
    Zwraca obiekt HydrogenBondAnalysis.
    """
    logger.info("Rozpoczynam analizę wiązań wodorowych...")
    # In newer versions of MDAnalysis, the parameters are positional rather than keyword
    h = HydrogenBondAnalysis(universe, 
                             "resname SOL or resname TIP4", # donor selection string
                             "resname SOL or resname TIP4", # acceptor selection string
                             d_h_a=dist, 
                             angle_dha=angle,
                             start=start, stop=stop, step=step)
    h.run()
    return h 