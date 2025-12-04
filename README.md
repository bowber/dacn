# Báo cáo Thí nghiệm Cơ sở Điều khiển

Báo cáo thí nghiệm môn học "Cơ sở điều khiển quá trình" - Khoa Kỹ thuật Hóa học, Trường Đại học Bách Khoa TP.HCM.

## Build

```bash
make pdf      # Build PDF using Docker (recommended)
make local    # Build locally (requires xelatex)
make clean    # Remove auxiliary files
```

## Cấu trúc

```
src/
├── main.tex      # Entry point
├── preamble.tex  # LaTeX preamble và styles
├── cover_page.tex
├── lab2.tex      # Bài TN 2: Thiết bị điều khiển và lập trình vi điều khiển
├── lab3.tex      # Bài TN 3: Khảo sát bộ điều khiển ON-OFF và PID
└── lab4.tex      # Bài TN 4: Thiết kế hệ thống điều chỉnh các đại lượng cơ bản
```

## Nội dung các bài thí nghiệm

### Bài thí nghiệm 2

Bài TN 2 tập trung vào **Thí nghiệm 3: Mô phỏng vận hành điều khiển hệ thống chưng cất**.

Nội dung chính:
- Điều khiển trình tự (Sequence Control) cho hệ thống chưng cất
- 3 chu trình: Khởi động (ON), Làm việc ổn định (SS), Dừng (OFF)
- Sử dụng Arduino để lập trình điều khiển

### Bài thí nghiệm 3

Khảo sát các bộ điều khiển ON-OFF và PID trên các hệ thống điều khiển đại lượng cơ bản.

### Bài thí nghiệm 4

Thiết kế bộ điều khiển PID theo phương pháp Ziegler-Nichols.

## Versioning

Versioned via git tags using semantic versioning: `vMAJOR.MINOR.PATCH`
