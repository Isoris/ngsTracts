# notes.md

Working notes and key insights for ngsTracts. Not normative — for normative
specs see `docs/METHODOLOGY.md` and `docs/SCHEMA.md`.

---

## Gene conversion segregation patterns (4:4 / 5:3 / 6:2)

These ratios are NOT about chromosome number. They are about how many DNA
strands / chromatids carry one allele vs the other allele after meiosis.

In fungal tetrad genetics, people count the 8 DNA products after meiosis +
one postmeiotic DNA replication. For a single heterozygous site:

| Pattern | Meaning |
| ------- | ------- |
| 4:4     | Normal Mendelian segregation — no gene conversion |
| 6:2     | Mismatch was repaired in favor of one allele — classic gene conversion signature |
| 5:3     | Mismatch was NOT repaired before postmeiotic segregation — unresolved heteroduplex |

### Why 6:2 happens

Imagine the homologs differ at one SNP:

- Blue homolog: `A`
- Red homolog:  `a`

During recombination, one strand from blue pairs with one strand from red.
This creates heteroduplex DNA — one strand says `A`, the other says `a`.
That mismatch is unstable.

If the repair system converts the red `a` into blue `A`, one chromatid
has been genetically converted. Instead of:

```
A A A A | a a a a    (4:4)
```

you get:

```
A A A A A A | a a    (6:2)
```

Six products carry `A`, two carry `a` — classic gene-conversion signature.

### Why 5:3 happens

5:3 happens when the mismatch is **not repaired** before the next DNA
replication / cell division. One DNA duplex remains mixed:

- one strand = `A`
- one strand = `a`

After replication, that mixed molecule gives one `A` product and one `a`
product. Final count:

```
A A A A A | a a a    (5:3)
```

This is called **postmeiotic segregation**.

### Link to SDSA / DSB repair

The DSB repair process creates a short region where one strand from
homolog A pairs with one strand from homolog B. If they differ by SNPs,
this region contains mismatches. Then:

- Mismatch repaired toward donor allele → gene conversion → **6:2**
- Mismatch left unrepaired → postmeiotic segregation → **5:3**

---

## Implications for the catfish WGS atlas

For our catfish WGS data, we usually will NOT literally observe 5:3 or 6:2,
because those ratios require tetrad-style observation or very controlled
family/meiotic-product assays. We have parent–offspring trios, not tetrads.

In our data, the equivalent signal is weaker: a short tract copied from one
arrangement into the other — a gene-conversion-like / SDSA-compatible
footprint.

### Labeling rule for the atlas

Do **not** label things "5:3" or "6:2" unless we have a real tetrad-style
assay. Use these labels instead:

- `conversion-like tract`
- `postmeiotic-segregation-like ambiguity`
- `unresolved heteroduplex-compatible pattern`
- `DSB-repair footprint`

These map onto the existing ngsTracts classes (see README, `docs/METHODOLOGY.md`):
the `NCO` class is the closest analogue — non-crossover gene conversion
inferred from a short departure interval with matching flanks.
