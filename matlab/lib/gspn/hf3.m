function key = hf3(marking)
    B = getByteStreamFromArray(marking);
    md = java.security.MessageDigest.getInstance('SHA-1');
    md.update(B);
    key = hex2dec(reshape(dec2hex(typecast(md.digest(),'uint8'))',1,[]));
end