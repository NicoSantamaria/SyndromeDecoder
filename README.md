# SyndromeDecoder

How can you decipher a message once errors have been introduced into it? In particular, if a 
message is transmitted through an imperfect channel — that is, a channel with a non-zero 
probability of introducing errors into the message — how can the original message be safely recovered?

Most well-known, effectives codes are linear codes. Linear codes are vector spaces over a finite 
field, where messages can be expressed as individual vectors. When an erroneous message is received, 
the goal is to decode the message to the "nearest" codeword within the vector space, with respect 
to a metric. Syndrome decoding is a technique that allows for messages to be decoded with minimal 
storage space required after initial computations, harnessing the fact that many vectors within a 
vector space belong to the the same cosets within that vector space, and those cosets can be represented 
by computationally inexpensive syndromes.

In Syndrome Decoder, the Python module MatrixMod.py handles the background mathematics (linear algebra 
operations over a finite field), while LinearCode.py performs the syndrome decoding algorithm. 
Syndrome Decoder is a Python program which uses syndrome decoding to correct for errors in linear codes.
