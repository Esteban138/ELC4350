
import binascii
if __name__ == '__main__':
    f = open("decoded_message.txt", "r")
    str = f.read()
    f.close()
    counter = 1;
    binary_int = int(str, 2)
    decoded_message =  binascii.unhexlify('%x' % binary_int)
    print(decoded_message)

    # USE FOR TESTING ONLY

    # for x in str:
    #     print(x, end='')
    #     if(counter %8 == 0 and counter !=0):
    #         print(" ", end='')
    #     counter = counter +1
