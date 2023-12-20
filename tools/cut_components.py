DEBUG=True
def cut_components(Nei,Cut,CutSize,Fal, isFalse):
  # Define the connected components of the Cut
  if DEBUG:
    print(f"Checking type and value of CutSize:")
    print(f"CutSize: {CutSize}")
    print(f"type(CutSize): {type(CutSize)}")
  if CutSize == [1]:
    # Cut is connected and therfore Study is also
    CompSize = [1]
    Components = {1: np.array(Cut).astype(int)}
    return Components, CompSize
  elif CutSize == [2]:
    I = Nei[Cut[0]] = Cut[1]
    if np.any(I):
      Components = {1:np.array(Cut).astype(int)}
      CompSize = [1]
    else:
      Components = {1:np.array(Cut[0]).astype(int), 2:np.array(Cut[1]).astype(int)}
      CompSize = [1, 1]
    return Components, CompSize
  elif CutSize == 3:
    I = Nei[Cut[0]] == Cut[1]
    J = Nei[Cut[0]] == Cut[2]
    K = Nei[Cut[1]] == Cut[2]
    if np.any(I) + np.any(J) + np.any(K) >=2:
      CompSize = [1]
      Components = {1:np.array(Cut).astype(int)}
    elif np.any(I):
      Components = {1:np.asarray(Cut[:2]).astype(int), 2:np.array(Cut[2]).astype(int)}
      CompSize = [2,1]
    elif np.any(J):
      Components = {1:np.asarray(Cut[0,2]).astype(int), 2:np.array(Cut[1]).astype(int)}
      CompSize = [2,1]
    elif np.any(K):
      Components = {1:np.asarray(Cut[1,2]).astype(int), 2:np.array(Cut[0]).astype(int)}
      CompSize = [2,1]
    else:
      CompSize[1, 1, 1]
      Components = {1:np.asarray(Cut[0]).astype(int), 2:np.array(Cut[1]).astype(int), 3:np.array(Cut[2]).astype(int)}
  else: 
    Components = {}
    CompSize = np.zeros(CutSize)
    Comp = np.zeros(CutSize)
    Fal[Cut] = True
    numC = 0
    m = Cut[0]
    i = 0
    while i < CutSize:
      Added = Nei[m]
      I = Fal[Added]
      Added = Added[I]
      a = len(Added)
      Comp[0] = m
      Fal[m] = False
      t = 1
      while a > 0:
        Comp[t+1-1:t+a] = Added
        Fal[Added] = False
        t = t + a
        Ext = np.concatenate([Nei[key] for key in Nei.keys()])
        Ext = unique_elemnts(Ext,isFalse)
        I = Fal[Ext]
        Added = Ext[I]
        a = len(Added)
      i += t
      numC +=1
      Components[numC] = Comp[:t]
      CompSize[numC] = t
      if i < CutSize:
        J = Fal[Cut]
        m = Cut[J]
        m = m[0]
    Components = Components[:numC]
    CompSize = CompSize[:numC]
    return Components, CompSize

