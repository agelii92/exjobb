import os
def main():
    lines = open('simulation_settings.txt','r').readlines()
    out = open('make.py','a')
    for i in lines:
        if not (len(i) == 1 or i[0] == '#'):
            out.write(i)
    with open('file.txt', 'a') as file:
        file.write('input')
    lines = open('code.py','r')
    for i in lines:
        out.write(i)

if __name__ == "__main__":
    main()
    exec(open("./make.py").read())
