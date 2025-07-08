# 🎯 Projeto Three.js - Visualização 3D

## 📋 Descrição do Projeto

Este é um projeto acadêmico de visualização 3D desenvolvido em JavaScript puro utilizando a biblioteca Three.js. O projeto demonstra conceitos fundamentais de computação gráfica 3D, incluindo renderização, shaders customizados, texturas, animações e diferentes tipos de câmera.

## 🛠️ Tecnologias Utilizadas

- **JavaScript ES6+** - Linguagem principal
- **Three.js (r128)** - Biblioteca para gráficos 3D WebGL
- **WebGL** - API de renderização 3D no navegador
- **GLSL** - Linguagem de shaders (vertex e fragment)
- **HTML5 Canvas** - Para geração de texturas procedurais
- **CSS3** - Estilização da interface

## 📁 Estrutura do Projeto

```
ProjetoPG/
├── index.html              # Página principal com canvas WebGL
├── main.js                 # Código principal da aplicação
├── shaders/                # Pasta com shaders customizados
│   ├── custom-vertex.glsl  # Vertex shader com animação de ondas
│   └── custom-fragment.glsl # Fragment shader com cores animadas
├── assets/                 # Recursos do projeto
│   └── texture-generator.js # Gerador de texturas procedurais
└── README.md              # Documentação do projeto
```

## ✅ Requisitos Implementados

### 🌟 Funcionalidades Básicas
- [x] **Cena 3D inicial** com fundo cinza escuro (`#1a1a1a`)
- [x] **Objetos 3D simples**: Cubo, esfera e plano
- [x] **Duas câmeras**: Perspectiva e ortográfica com alternância
- [x] **Shader customizado** aplicado ao cubo usando `RawShaderMaterial`
- [x] **Texturas básicas** aplicadas em objetos (gradiente e xadrez)
- [x] **Animações** de rotação e translação nos objetos
- [x] **Posicionamento e escala** individuais dos objetos
- [x] **Controles de câmera** com mouse usando `OrbitControls`

### 🎨 Detalhes da Implementação

#### 1. **Cena e Iluminação**
- Fundo cinza escuro para contraste
- Luz ambiente suave (0x404040)
- Luz direcional principal (0xffffff)

#### 2. **Objetos 3D**
- **Cubo**: Material com shader customizado, posição (-3, 0, 0), escala 1.2x
- **Esfera**: Textura gradiente radial, posição (3, 0, 0), escala 1.0x, movimento de quique
- **Plano**: Textura xadrez 32x32, posição (0, -2, 0), rotacionado horizontalmente

#### 3. **Câmeras**
- **Perspectiva**: FOV 75°, posição (5, 5, 5)
- **Ortográfica**: Frustum 10x10, mesma posição
- Alternância com tecla `C`

#### 4. **Shader Customizado**
- **Vertex Shader**: Cria ondulações baseadas em seno/cosseno
- **Fragment Shader**: Cores animadas baseadas na posição e tempo
- Uniforms: `time`, `projectionMatrix`, `modelViewMatrix`

#### 5. **Texturas Procedurais**
- **Gradiente**: Cores do amarelo ao rosa (#ffeb3b → #e91e63)
- **Xadrez**: Padrão alternado vermelho/turquesa (#ff6b6b / #4ecdc4)

## 🎮 Como Interagir com a Cena

### 🖱️ Controles do Mouse
- **Botão esquerdo + arrastar**: Rotacionar câmera ao redor da cena
- **Botão direito + arrastar**: Mover câmera (pan)
- **Scroll**: Zoom in/out
- **Botão do meio + arrastar**: Zoom dinâmico

### ⌨️ Controles do Teclado
- **Tecla C**: Alternar entre câmera perspectiva e ortográfica

### 📱 Responsividade
- A aplicação se adapta automaticamente ao redimensionamento da janela
- Ambas as câmeras mantêm proporção correta

## 🎭 Animações Implementadas

1. **Cubo com Shader**:
   - Rotação contínua nos eixos X e Y
   - Ondulações procedurais na superfície
   - Cores que mudam dinamicamente

2. **Esfera com Textura**:
   - Movimento de quique realista (não atravessa o plano)
   - Rotação no eixo Y
   - Textura gradiente estática

3. **Plano Base**:
   - Rotação suave no eixo Z
   - Textura xadrez com padrão repetido
   - Tamanho expandido (32x32 unidades)

## 🚀 Como Executar

1. **Clone ou baixe o projeto**
2. **Abra o arquivo `index.html`** em um navegador moderno
3. **Aguarde o carregamento** dos recursos Three.js via CDN
4. **Interaja com a cena** usando mouse e teclado

### 📋 Requisitos do Sistema
- Navegador com suporte WebGL (Chrome, Firefox, Safari, Edge)
- Conexão com internet (para CDN do Three.js)
- Hardware com suporte básico para aceleração gráfica

## 🔧 Arquitetura do Código

### 📦 Módulos Principais

1. **Inicialização** (`init()`):
   - Carregamento assíncrono de shaders
   - Configuração da cena, câmeras e renderer
   - Criação de objetos e configuração de eventos

2. **Gerenciamento de Shaders** (`loadShaders()`):
   - Carregamento dinâmico de arquivos GLSL
   - Fallback para shaders inline em caso de erro

3. **Criação de Objetos** (`createObjects()`):
   - Factory methods para cada tipo de objeto
   - Aplicação de materiais e texturas específicos

4. **Loop de Animação** (`animate()`):
   - Renderização contínua a 60 FPS
   - Atualização de animações baseadas no tempo
   - Sincronização com `requestAnimationFrame`

## 📊 Performance e Otimizações

- **Damping nos controles**: Suaviza movimentos da câmera
- **Reutilização de texturas**: Canvas gerados uma única vez
- **Uniforms otimizados**: Apenas `time` é atualizado por frame
- **Geometrias simples**: Baixo número de vértices para performance

## 🎓 Valor Acadêmico

Este projeto demonstra:
- **Fundamentos de WebGL** através do Three.js
- **Programação de shaders** em GLSL
- **Transformações 3D** (rotação, translação, escala)
- **Sistemas de coordenadas** e projeções
- **Texturas procedurais** e mapeamento UV
- **Iluminação básica** (ambiente + direcional)
- **Controles de câmera** e interatividade
- **Arquitetura modular** de aplicações gráficas

## 🐛 Troubleshooting

### Problema: Shaders não carregam
**Solução**: O código possui fallback automático para shaders inline

### Problema: Texturas não aparecem
**Solução**: Verifique se o arquivo `texture-generator.js` está carregando

### Problema: Controles não funcionam
**Solução**: Certifique-se de que o OrbitControls está carregando via CDN

### Problema: Performance baixa
**Solução**: Teste em outro navegador ou verifique suporte WebGL

---

**Desenvolvido para fins acadêmicos** | **Three.js r128** | **WebGL** 