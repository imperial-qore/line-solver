def allbut(y = None,xset = None):
    # Y=ALLBUT(Y,XSET)

    y = y(difference(arange(0,len(y)+1),xset))
    return y