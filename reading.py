def read(filename):
    with open(filename, 'r') as file:
        dados = file.readlines()
    dados = [line.strip().split() for line in dados]
    count=0
    nnodes=int(dados[count][1])
    count+=1
    nmats=int(dados[count][1])
    count+=1
    nthick=int(dados[count][1])
    count+=1
    nelem3=int(dados[count][1])
    count+=1
    nelem6=int(dados[count][1])
    count+=1
    nelem10=int(dados[count][1])

    nelems=nelem3+nelem6+nelem10

    count+=1
    nforces=int(dados[count][1])
    count+=1
    npressures=int(dados[count][1])
    count+=1
    ndisplas=int(dados[count][1])
    count+=3

    return dados,nnodes,nmats,nthick, nelem3, nelem6, nelem10, nelems,nforces,npressures,ndisplas,count
