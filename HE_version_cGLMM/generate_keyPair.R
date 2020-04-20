library(homomorpheR)

keyPair <- PaillierKeyPair$new(modulusBits = 1024)

pub_key=keyPair$pubkey
priv_key = keyPair$getPrivateKey()

print(pub_key)
print(priv_key)
save(pub_key, file="pub_key.RData")
save(priv_key, file="priv_key.RData")


