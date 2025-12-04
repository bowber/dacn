# Báo cáo Thí nghiệm Thiết bị Đo lường và Điều khiển

Báo cáo thí nghiệm môn học "Thiết bị Đo lường và Điều khiển" (CH4132) - Khoa Kỹ thuật Hóa học, Trường Đại học Bách Khoa, ĐHQG-HCM.

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
├── lab1.tex       # Bài TN 1: Khảo sát cảm biến nhiệt độ
├── lab2.tex       # Bài TN 2: Điều khiển ON/OFF
├── lab3.tex       # Bài TN 3: Điều khiển PID
└── lab4.tex       # Bài TN 4: Factory I/O

data/              # Dữ liệu thí nghiệm
assets/            # Hình ảnh, logo
docs/              # Tài liệu tham khảo
output/            # PDF output
```

## Thiết bị và phần mềm

- **PLC:** Siemens S7-1200
- **Phần mềm lập trình:** TIA Portal V16
- **Mô phỏng:** Factory I/O

## Nội dung các bài thí nghiệm

### Bài thí nghiệm 1: Khảo sát cảm biến nhiệt độ

Khảo sát đặc tính và hiệu chuẩn cảm biến nhiệt độ.

- Cảm biến nhiệt điện trở (RTD, Pt100)
- Cặp nhiệt điện (Thermocouple)
- Đường cong hiệu chuẩn

### Bài thí nghiệm 2: Điều khiển ON/OFF

Thiết kế và khảo sát hệ thống điều khiển nhiệt độ ON/OFF sử dụng PLC S7-1200.

- Lập trình PLC với TIA Portal
- Điều khiển bang-bang
- Phân tích dao động và hysteresis

### Bài thí nghiệm 3: Điều khiển PID

Thiết kế và chỉnh định bộ điều khiển PID cho hệ thống điều khiển nhiệt độ.

- Cấu hình khối PID trong TIA Portal
- Chỉnh định thông số Kp, Ki, Kd
- Đánh giá chất lượng điều khiển

### Bài thí nghiệm 4: Factory I/O

Mô phỏng và điều khiển hệ thống công nghiệp với Factory I/O.

- Kết nối PLC S7-1200 với Factory I/O
- Lập trình điều khiển quy trình tự động
- Giám sát và vận hành hệ thống ảo

## Versioning

Versioned via git tags using semantic versioning: `vMAJOR.MINOR.PATCH`
