#!/usr/bin/env julia
#
# Verifies that every tracked `.jl`, `.md` and `.toml` file in this repository is Unicode
# NFC-normalised.
#
#     julia --startup-file=no scripts/verify_nfc.jl
#
# Exit status is 0 when every file is normalised and 1 otherwise, so this is usable as a gate.
#
# Why the invariant is worth asserting. NFD normalisation leaves the letters
# `ū`, `ḡ`, `ṗ` and `ẋ` stored as a base letter plus a combining mark rather than as one
# codepoint. That makes no difference to the compiled code -- Julia's parser normalises
# identifiers to NFC, so an NFD source file produces byte-identical symbols -- but it defeats
# every byte-matching tool. A `grep`, an editor search or an automated replacement typed in NFC
# matches nothing in an NFD file, silently.
#
# String literals are the exception that makes the invariant worth keeping rather than merely
# tidy: they are not parser-normalised, so in an NFD file `:ẋ` and `Symbol("ẋ")` compare unequal
# while looking identical. The `Base.show` methods on the DAE types print literal `"   ū = "` and
# `"   ḡ = "`, so their emitted bytes depend on this invariant holding.
#
# Note that NFC does not eliminate every combining mark. `q̇`, `v̄` and `f̄` have no precomposed
# codepoint and remain two codepoints; the assertion below is `s == normalize(s, :NFC)`, not the
# absence of combining marks.

using Unicode

cd(dirname(@__DIR__))

const FILES = filter(!isempty, split(readchomp(`git ls-files -- "*.jl" "*.md" "*.toml"`), '\n'))

failures = String[]

for path in FILES
    isfile(path) || continue
    source = read(path, String)
    if !isvalid(source)
        push!(failures, "$path: not valid UTF-8")
        continue
    end
    source == Unicode.normalize(source, :NFC) || push!(failures, "$path: not NFC")
end

println(length(FILES), " tracked source files checked")

if isempty(failures)
    println("all NFC-normalised")
    exit(0)
else
    println(length(failures), " failure(s):")
    foreach(f -> println("    ", f), failures)
    exit(1)
end
