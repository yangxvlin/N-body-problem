if __name__ == "__main__":
    datas = [(8, 90), (50, 90), (100, 90)]
    for n, k in datas:
        with open("nqueen_{}_{}.data".format(n, k), 'w') as f:
            print(n, file=f)
            print(k, file=f)
