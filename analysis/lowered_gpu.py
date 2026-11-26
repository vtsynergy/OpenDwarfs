"""Static lowered-code features of GPU kernels from a built binary (cuobjdump): resource usage + SASS instruction-class histogram."""
import re, subprocess, collections, shutil, os
CUOBJDUMP = shutil.which("cuobjdump") or "/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/compilers/bin/cuobjdump"
CLASSES = [("mem_global", r"^(LDG|STG|ATOMG|RED|LD|ST|ATOM)\b"), ("mem_shared", r"^(LDS|STS|ATOMS|LDSM)\b"), ("mem_local", r"^(LDL|STL)\b"),
           ("mem_const", r"^(LDC|ULDC)\b"), ("fp", r"^(FADD|FMUL|FFMA|FMNMX|FSETP|FSEL|DADD|DMUL|DFMA|DSETP|MUFU|HADD2|HFMA2|FRND|F2F|F2I|I2F|I2FP)\b"),
           ("sync", r"^(BAR|WARPSYNC|BSYNC|BSSY|MEMBAR|ERRBAR|DEPBAR|YIELD|NANOSLEEP|VOTE|SHFL)\b"), ("control", r"^(BRA|BRX|JMP|CALL|RET|EXIT|BREAK|KILL|JMX|BPT)\b"),
           ("int", r"^(IADD3|IADD|IMAD|IMUL|LOP3|LOP|SHF|SHL|SHR|ISETP|LEA|SEL|MOV|UMOV|IABS|PRMT|POPC|FLO|BFE|BFI|ISCADD|IMNMX|UIADD3|ULOP3|USHF|UISETP|ULEA|S2R|S2UR|CS2R|R2UR|P2R|R2P|PLOP3|UPLOP3|NOP)\b")]
def _cls(op):
    for n, rx in CLASSES:
        if re.match(rx, op): return n
    return "other"
def kernel_key(name):
    """Normalise clang (__omp_offloading_.._func_lNN) / nvhpc (nvkernel_func_FnLnn_k) / Itanium-mangled kernel names to the source function name."""
    m = re.match(r"__omp_offloading_[0-9a-f]+_[0-9a-f]+_(.*)_l\d+$", name)
    inner = m.group(1) if m else name
    m = re.match(r"nvkernel_(.*)_F\d+L\d+_\d+$", name)
    if m: inner = m.group(1)
    if inner.startswith("_Z"):
        s = inner[2:]
        if s[:1] in "LN": s = s[1:]
        last = None
        while re.match(r"\d+", s):
            n = int(re.match(r"\d+", s).group(0)); s = s[len(str(n)):]; last = s[:n]; s = s[n:]
        if last: return last
    return inner
def extract(binary):
    """-> {key: dict(regs, shared, local, const, instr=int, classes={cls:n}, raw_name)}; {} if the binary has no device code."""
    out = {}
    try:
        res = subprocess.run([CUOBJDUMP, "-res-usage", binary], capture_output=True, text=True, timeout=120).stdout
        sass = subprocess.run([CUOBJDUMP, "-sass", binary], capture_output=True, text=True, timeout=120).stdout
    except Exception:
        return out
    for m in re.finditer(r"Function (\S+):\s*\n\s*(REG:\d+.*)", res):
        d = dict(re.findall(r"(\w+(?:\[\d+\])?):(\d+)", m.group(2)))
        out[kernel_key(m.group(1))] = dict(raw_name=m.group(1), regs=int(d.get("REG", 0)), stack=int(d.get("STACK", 0)), shared=int(d.get("SHARED", 0)),
                                          local=int(d.get("LOCAL", 0)), instr=0, classes={})
    cur = None
    for line in sass.splitlines():
        m = re.match(r"\s*Function : (\S+)", line)
        if m: cur = kernel_key(m.group(1)); out.setdefault(cur, dict(raw_name=m.group(1), regs=0, stack=0, shared=0, local=0, instr=0, classes={})); continue
        m = re.match(r"\s*/\*[0-9a-f]{4}\*/\s+(?:@!?U?P\w\s+)?([A-Z][A-Z0-9_]*)", line)
        if m and cur:
            c = _cls(m.group(1)); k = out[cur]; k["instr"] += 1; k["classes"][c] = k["classes"].get(c, 0) + 1
    return out
