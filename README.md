# Đồ án Chuyên ngành - Tối ưu hóa năng suất sản phẩm đỉnh tháp chưng cất ethanol-nước

Đồ án Chuyên ngành ngành Kỹ thuật Hóa học - Khoa Kỹ thuật Hóa học, Trường Đại học Bách Khoa, ĐHQG-HCM.

**Đề tài:** Tối ưu hóa năng suất sản phẩm đỉnh của tháp chưng cất ethanol-nước thông qua điều khiển tự động áp suất đỉnh tháp.

## Tóm tắt

Đồ án nghiên cứu và triển khai hệ thống điều khiển áp suất đỉnh tháp chưng cất ethanol-nước quy mô phòng thí nghiệm. Thông qua việc điều khiển lưu lượng nước làm mát bằng bộ điều khiển PID, hệ thống có thể:

- Giảm lượng nước làm mát từ 7.2 L/min xuống 3.8 L/min (tiết kiệm 47%)
- Tăng năng suất sản phẩm đỉnh nhờ giảm ngưng tụ quá mức
- Duy trì chất lượng sản phẩm (90% vol ethanol)

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
├── main.tex          # Entry point
├── preamble.tex      # LaTeX preamble và styles
├── cover_page.tex    # Trang bìa
├── assignment.tex    # Nhiệm vụ đồ án
├── acknowledgements.tex # Lời cảm ơn
├── abstract.tex      # Tóm tắt
├── chapter1.tex      # Chương 1: Tổng quan
├── chapter2.tex      # Chương 2: Cơ sở lý thuyết
├── chapter3.tex      # Chương 3: Mô hình thí nghiệm
├── chapter4.tex      # Chương 4: Tính toán và thiết kế
├── chapter5.tex      # Chương 5: Kết luận và kế hoạch
└── references.tex    # Tài liệu tham khảo

assets/               # Hình ảnh thiết bị
docs/                 # Tài liệu tham khảo
simulation/           # Mô phỏng Python
output/               # PDF output
```

## Thiết bị chính

| Thiết bị | Model | Chức năng |
|----------|-------|-----------|
| PLC | Siemens S7-1200 | Điều khiển trung tâm |
| Cảm biến áp suất | SITRANS P200 | Đo áp suất đỉnh tháp |
| Cảm biến lưu lượng | Kobold DPM | Đo lưu lượng hồi lưu |
| Cảm biến mức | Kobold NMC | Đo mức bình hồi lưu |
| Van điều khiển | Bürkert Type 2873 + 8605 | Điều khiển nước làm mát |
| Reboiler | Điện trở 2×3kW | Gia nhiệt đáy tháp |

## Thông số hệ thống

- Tháp chưng cất: 6 mâm xuyên lỗ, đường kính 150mm
- Công suất reboiler: 6000W (2×3000W)
- Condenser: Ống xoắn ruột gà, diện tích 0.28 m²
- Nồng độ nhập liệu: 10% vol
- Nồng độ sản phẩm đỉnh: 90% vol

## Versioning

Versioned via git tags using semantic versioning: `vMAJOR.MINOR.PATCH`
