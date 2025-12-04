# Báo cáo Thí nghiệm Cơ sở Điều khiển Quá trình

Báo cáo thí nghiệm môn học "Cơ sở điều khiển quá trình" (MSMH 3342) - Khoa Kỹ thuật Hóa học, Trường Đại học Bách Khoa, ĐHQG-HCM.

## Build

```bash
make pdf      # Build PDF using Docker (recommended)
make local    # Build locally (requires xelatex)
make clean    # Remove auxiliary files
make distclean # Remove all generated files including PDFs
make shell    # Open interactive shell in Docker for debugging
```

## Cấu trúc dự án

```
src/
├── main.tex       # Entry point
├── preamble.tex   # LaTeX preamble và styles
├── cover_page.tex # Trang bìa
├── lab2.tex       # Bài TN 2: Thiết bị điều khiển và lập trình vi điều khiển
├── lab3.tex       # Bài TN 3: Khảo sát bộ điều khiển ON-OFF và PID
└── lab4.tex       # Bài TN 4: Khảo sát ảnh hưởng các thông số bộ điều khiển PID

data/              # Dữ liệu thí nghiệm (.csv, .TRD)
assets/            # Hình ảnh, logo
docs/              # Tài liệu tham khảo (PDF)
output/            # PDF output
```

## Nội dung các bài thí nghiệm

### Bài thí nghiệm 2: Thiết bị điều khiển và lập trình vi điều khiển

Mô phỏng vận hành điều khiển hệ thống chưng cất sử dụng Arduino.

- Điều khiển trình tự (Sequence Control)
- 3 chu trình: Khởi động (ON), Làm việc ổn định (SS), Dừng (OFF)
- Lập trình Arduino điều khiển bơm, van, điện trở

### Bài thí nghiệm 3: Khảo sát bộ điều khiển ON-OFF và PID

Khảo sát và so sánh các bộ điều khiển trên hệ thống điều khiển mức chất lỏng với phần mềm Halalab.

- Bộ điều khiển ON-OFF
- Bộ điều khiển P, PI, PID
- Đánh giá chỉ tiêu chất lượng: thời gian xác lập, độ quá điều chỉnh, sai lệch tĩnh

### Bài thí nghiệm 4: Khảo sát ảnh hưởng các thông số bộ điều khiển PID

Khảo sát ảnh hưởng của các thông số $K_P$, $K_I$, $K_D$ đến chất lượng điều khiển.

- Thay đổi từng thông số và ghi nhận đáp ứng
- Xác định bộ thông số tối ưu

## Tài liệu tham khảo

1. B. N. Pha, *Điều khiển Quá trình Công nghệ Hoá học -- Cơ sở Điều khiển Quá trình*, Quyển 1. Trường Đại học Bách khoa, ĐHQG-HCM, 2021.
2. B. N. Pha và P. H. H. P. Lợi, *Điều khiển Quá trình Công nghệ Hoá học -- Hướng dẫn Thí nghiệm – Thực hành Cơ sở Điều khiển Quá trình*, Quyển 2. Trường Đại học Bách khoa, ĐHQG-HCM, 2021.

## Versioning

Versioned via git tags using semantic versioning: `vMAJOR.MINOR.PATCH`
