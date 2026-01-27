from cli import parse_args
from build_config import build_config
from validate_config import validate_config
from pathlib import Path
import sys, yaml

def main():
    args = parse_args()
    
    try:
        cfg = build_config(args)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    errors = validate_config(cfg)
    if errors:
        print("Refusing to write invalid config:")
        for e in errors:
            print(f" - {e}", file=sys.stderr)
        sys.exit(1)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    yaml.safe_dump(cfg, out.open("w"), sort_keys=False)
    print(f"Wrote config to {out}")

if __name__ == "__main__":
    main()

